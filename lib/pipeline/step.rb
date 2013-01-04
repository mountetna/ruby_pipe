require 'fileutils'

module Pipeline
  module Step
    include Pipeline::Logger
    # This creates a new step in a pipeline
    module ClassMethods
      def runs_tasks(*tasklist)
        @tasks = tasklist
      end
      def tasks
        @tasks
      end
      def required
        @required ||= @tasks.map{ |t| self.class.daughter_class(t).required_files }.flatten
      end
      def made
        @made ||= @tasks.map{ |t| self.class.daughter_class(t).made_files }.flatten
      end
      def input
        required - made
      end
    end
    def self.included(base)
      base.extend(ClassMethods)
    end

    def initialize(s)
      @script = s
      config.set_opt :step, step_name
    end

    def step_name; self.class.name.split(/::/).last.snake_case.to_sym; end
    def script; @script; end
    def config; @script.config; end

    def schedule_job(action,opts=nil)
      # there are some standard vars to pass in
      vars = {
        :CONFIG => config.config_file,
        :STEP => config.step,
        :ACTION => action,
        :LIB_DIR => config.lib_dir,
        :SINGLE_STEP => config.single_step
      }
      opts = { 
        :N => "#{config.pipe}.#{config.step}",
        :m => "n",
        :j => "oe",
        :o => "test.log",
        :l => "walltime=3:00:00:00",
        :v => vars.map{|v,n| "#{v}=#{n}" }.join(",")
      }.merge(opts||{})


      `/opt/moab/bin/msub #{opts.map{ |o,v| "-#{o} #{v}" }.join(" ")} #{config.lib_dir}/run_step.rb |tail -1`
    end

    def setup_scheduler(prevjob)
      # setup the scheduler to run this task.
      setup_logging

      log_info "Scheduling #{step_name}"
      schedule_job(:schedule)
    end

    def setup_exec
      # setup the scheduler to execute this task.
      setup_logging

      log_info "Starting execution for #{step_name}".yellow.bold
      if config.splits
        schedule_job(:exec, :t => "1-#{config.splits}")
      else
        schedule_job(:exec)
      end
    end

    def complete
      # consolidate logs for this step. 
      
      # Error out if there is an error pid
      if File.exists?(config.error_pid)
        FileUtils.rm(config.error_pid)
        File.open(config.error_file, "w") do |f|
          f.puts "Script failed at #{config.prev_step}, see logs in #{config.log_dir}/. Resume with 'run_pipeline start #{config.config_file} #{config.prev_step}'"
        end
      end
    end

    def cleanup
      return if config.keep_temp_files
      self.class.tasks.each do |t|
        self.class.daughter_class(t).dump_files.each do |f|
          filename = config.send(f)
          if filename && filename.is_a?(Array)
            filename.each do |fn|
              FileUtils.rm(fn) if fn && File.size?(fn) && File.writable?(fn)
            end
          else
            FileUtils.rm(filename) if filename && File.size?(filename) && File.writable?(filename)
          end
        end
      end
    end

    def exec
      setup_logging
      self.class.tasks.each do |t|
        task = create_task t

        task.exec if task.should_run
      end
      return true
    end

    def audit
      log_info "Auditing #{self.class.name.snake_case} trial #{config.job_index}".yellow.bold

      self.class.tasks.each do |t|
        create_task(t).audit
      end
    end

    def create_task(t)
      self.class.daughter_class(t.to_sym).new(self)
    end
  end
end
