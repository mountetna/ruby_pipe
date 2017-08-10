require 'open3'
module Pipeline
  module Script
    def scheduler_type
      case
      when !`pgrep moab`.empty?
        :moab_scheduler
      when !`pgrep maui`.empty?
        :maui_scheduler
      when `qstat -help` =~ /^SGE/
        :sge_scheduler
      end
    end
  end
  module Scheduling
    class TorqueScheduler
      def run_job vars, opts
        opts = {
          :nodes => "1",
          :threads => "1",
          :memory => "1gb"
        }.merge(opts)

        
        fields = {
          :N => opts[:name],
          :W => opts[:wait] ? wait_string(opts[:wait],opts[:prev_trials]) : nil,
          :t => opts[:trials] ? "1-#{opts[:trials]}" : nil,
          :m => opts[:email_opts] ? opts[:email_opts] : "n",
          :M => opts[:email_addr] ? opts[:email_addr] : nil,
          :j => "oe",
          :o => "log/uncaught_errors.log",
          :v => vars.map{|v,n| n ? "#{v}=#{n}" : nil }.compact.join(",")
        }
        fields.delete_if { |k,v| v.nil? }
        

        res = []
        res.push "walltime=#{opts[:walltime]}:00:00:00" if opts[:walltime]
        res.push "nodes=#{opts[:nodes]}:ppn=#{opts[:threads]}"
        res.push "vmem=#{opts[:memory]}"

        yield res, fields
      end
    end

    class MoabScheduler < TorqueScheduler
      def wait_string job, trials
        "x=depend:afterany:#{job}"
      end

      def run vars, opts
        run_job vars,opts do |res,fields|
          cmd = "/opt/moab/bin/msub #{res.map{|r| "-l #{r}"}.join(" ")} #{fields.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{vars[:LIB_DIR]}/step_pipe.rb |tail -1"
          yield cmd if block_given?
          %x{#{cmd}}
        end
      end

      def cancel_id id
        system "canceljob '#{id}'"
      end

      def cancel(conf, allow_completion=nil)
        q = Nokogiri::XML(`showq --format=xml`)

        # in this case we want to find all of the matching jobs.
        q.xpath("//job[@User='#{ENV['USER']}']").each do |job|
          jid = job["JobID"].gsub(/\(.*\)/,"")
          # do more research on this JobID
          jobinfo = `checkjob -v -v "#{jid}"`

          if jobinfo =~ /Job Array Info/
            jobinfo = `checkjob -v -v "#{jid}[1]"`
          end

          jobinfo.scan(/(?:#PBS -v |EnvVariables:\s+)(.*)\n/) do |v|
            vars = Hash[v.first.split(/,/).map{|s| s.split(/=/)}]
            next if vars["CONFIG"] != conf
            next if vars["ACTION"] == "exec" && allow_completion
            # okay, you have a variable and a config file, cancel this job.
            cancel_id jid
          end
        end
      end
    end

    class MauiScheduler < TorqueScheduler
      def wait_string job, trials
        if job =~ /\[\]/
          "depend=afteranyarray:#{job}"
        else
          "depend=afterany:#{job}"
        end
      end

      def run vars, opts
        run_job vars,opts do |res,fields|
          fields[:o] = "/dev/null"
          cmd = "qsub #{res.map{|r| "-l #{r}"}.join(" ")} #{fields.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{vars[:LIB_DIR]}/step_pipe.rb |tail -1"
          %x{#{cmd}}
        end
      end

      def cancel_id ids
        system "qdel #{ids.join(" ")}"
      end

      def cancel(conf, allow_completion=nil)
        q = Nokogiri::XML(`qstat -f -x`)
        schedule_id = q.xpath("//Job[contains(Variable_List,'CONFIG=#{conf}') and contains(Variable_List,'ACTION=schedule')]/Job_Id").map(&:text)
        exec_id = q.xpath("//Job[contains(Variable_List,'CONFIG=#{conf}') and contains(Variable_List,'ACTION=exec')]/Job_Id").map(&:text)

        cancel_id (allow_completion ? schedule_id : schedule_id + exec_id)
      end
    end

    class GridScheduler
      def run_job vars, opts
        fields = {
          :N => opts[:name],
          :t => opts[:trials] ? "1-#{opts[:trials]}" : nil,
          :hold_jid => opts[:wait],
          :j => "y",
          :pe => opts[:threads] ? "alloc #{opts[:threads]}" : nil,
          :o => "log/uncaught_errors.log",
          :q => opts[:queue],
          :v => vars.map{|v,n| n ? "#{v}=#{n}" : nil }.compact.join(",")
        }
        fields.delete_if { |k,v| v.nil? }

        res = []
        res.push "internet=1" if opts[:internet]

        yield res, fields
      end
    end

    class SgeScheduler < GridScheduler
      def run vars, opts
        run_job vars,opts do |res,fields|
          fields[:o] = "/dev/null"
          ruby = %x{ readlink -f $(which ruby) }.chomp
          cmd = "/common/sge/bin/lx24-amd64/qsub -terse #{res.map{|r| "-l #{r}"}.join(" ")} #{fields.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{vars[:LIB_DIR]}/wrapper.sh #{vars[:LIB_DIR]}/step_pipe.rb"

          job = nil
          Open3.popen3(cmd) do |sin,sout,serr,wait_thr|
            job = sout.read.sub(/\..*/,'')
          end
          job.strip
        end
      end

      def cancel_id ids
        system "qdel #{ids.join(" ")}"
      end

      def cancel(conf, allow_completion=nil)
        q = Nokogiri::XML(`qstat -f -x`)
        schedule_id = q.xpath("//Job[contains(Variable_List,'CONFIG=#{conf}') and contains(Variable_List,'ACTION=schedule')]/Job_Id").map(&:text)
        exec_id = q.xpath("//Job[contains(Variable_List,'CONFIG=#{conf}') and contains(Variable_List,'ACTION=exec')]/Job_Id").map(&:text)

        cancel_id (allow_completion ? schedule_id : schedule_id + exec_id)
      end
    end

    def scheduler
      @scheduler ||= Pipeline::Scheduling.const_get(config.scheduler.to_s.camel_case).new
    end

    def cancel_job(allow_completion=nil)
      # find things matching the current job.
      scheduler.cancel config.config_file, allow_completion
    end

    def schedule_job(action,opts=nil)
      # there are some standard vars to pass in

      vars = {
        :CONFIG => config.config_file,
        :STEP => config.step,
        :PIPE => config.pipe,
        :SCRIPT => config.script,
        :ACTION => action,
        :LIB_DIR => config.lib_dir,
        :SCHEDULER => config.scheduler,
        :SINGLE_STEP => config.single_step
      }

      opts = {
        :name => "#{config.pipe}.#{action == :schedule ? :schedule : config.step}",
        :queue => config.step_queue,
        :email_addr => config.email_addr,
        :email_opts => config.email_opts
      }.merge(opts || {})
      
      scheduler.run(vars, opts) do |cmd|
        log_info "Scheduled #{config.step} with #{cmd}"
      end
    end
  end
end
