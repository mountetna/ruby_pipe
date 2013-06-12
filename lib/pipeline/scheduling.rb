module Pipeline
  module Script
    def scheduler_type
      case
      when !`pgrep moab`.empty?
        :moab_scheduler
      when !`pgrep maui`.empty?
        :maui_scheduler
      end
    end
  end
  module Scheduling
    class TorqueScheduler
      def run_job vars, opts
        fields = {
          :N => opts[:name],
          :W => opts[:wait] ? wait_string(opts[:wait],opts[:prev_splits]) : nil,
          :t => opts[:splits] ? "1-#{opts[:splits]}" : nil,
          :m => "n",
          :j => "oe",
          :o => "log/uncaught_errors.log",
          :v => vars.map{|v,n| n ? "#{v}=#{n}" : nil }.compact.join(",")
        }
        fields.delete_if { |k,v| v.nil? }

        res = []
        res.push "walltime=#{opts[:walltime]}:00:00:00" if opts[:walltime]
        res.push "nodes=1:ppn=#{opts[:threads]}" if opts[:threads]
        res.push "pmem=1gb"

        yield res, fields
      end
    end

    class MoabScheduler < TorqueScheduler
      def wait_string job, splits
        "x=depend:afterany:#{job}"
      end

      def run vars,opts
        run_job vars,opts do |res,fields|
        `/opt/moab/bin/msub #{res.map{|r| "-l #{r}"}.join(" ")} #{fields.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{vars[:LIB_DIR]}/step_pipe.rb |tail -1`
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
      def wait_string job, splits
        if job =~ /\[\]/
          "depend=afteranyarray:#{job}"
        else
          "depend=afterany:#{job}"
        end
      end

      def run vars, opts
        run_job vars,opts do |res,fields|
          fields[:o] = "/dev/null"
          `qsub #{res.map{|r| "-l #{r}"}.join(" ")} #{fields.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{vars[:LIB_DIR]}/step_pipe.rb |tail -1`
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
      }.merge(opts || {})

      scheduler.run(vars, opts)
    end
  end
end
