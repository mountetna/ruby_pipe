module Pipeline
  module Scheduling
    class MoabScheduler
      def run vars, opts
      fields = {
        :N => opts[:name],
        :W => opts[:wait] ? "x=depend:afterany:#{opts[:wait]}" : nil,
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

      `/opt/moab/bin/msub #{res.map{|r| "-l #{r}"}.join(" ")} #{fields.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{vars[:LIB_DIR]}/step_pipe.rb |tail -1`
      end
    end

    def scheduler_type
      case
      when !`pgrep moab`.empty?
        :moab_scheduler
      end
    end

    def scheduler
      @scheduler ||= Pipeline::Scheduling.const_get(config.scheduler.camel_case).new
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
        :SCHEDULER => config.scheduler
        :SINGLE_STEP => config.single_step
      }

      opts = { 
        :name => "#{config.pipe}.#{action == :schedule ? :schedule : config.step}",
        :walltime => 3,
        :threads => resources[:threads]
      }.merge(opts || {})

      scheduler.run(vars, opts)
    end
  end
end
