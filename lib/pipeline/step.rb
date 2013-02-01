require 'fileutils'

module Pipeline
  module Step
    include Pipeline::Logger
    include Pipeline::Scheduling
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

      def resources(opts = nil)
        @resources ||= {}
        @resources.update(opts || {})
      end
    end
    def self.included(base)
      base.extend(ClassMethods)
    end

    def initialize(s)
      @script = s
      config.set_opt :step, step_name
      setup_logging unless config.action == :audit
    end

    def step_name; self.class.name.split(/::/).last.snake_case.to_sym; end
    def script; @script; end
    def config; @script.config; end
    def resources; self.class.resources; end

    def setup_scheduler(prevjob)
      # setup the scheduler to run this task.
      prevjob.strip!

      log_info "Scheduling #{step_name}"
      schedule_job :schedule, :wait => prevjob
    end

    def setup_exec
      # setup the scheduler to execute this task.

      log_info "Starting execution for #{step_name}".yellow.bold
      if config.splits
        schedule_job :exec, :splits => config.splits
      else
        schedule_job :exec
      end
    end

    def complete
      # consolidate logs for this step. 
      
      # Error out if there is an error pid
      if File.exists?(config.error_pid)
        FileUtils.rm(config.error_pid)
        File.open(config.error_file, "w") do |f|
          f.puts "Script failed at #{config.step}, see logs in #{config.log_dir}/. Resume with '#{config.pipe}_#{config.script} start #{config.config_file} #{config.step}'"
        end
        exit
      end
    end

    def vacuum
    end

    def cleanup
      # run the vacuum script
      vacuum
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
      self.class.tasks.each do |t|
        task = create_task t

        task.exec if task.should_run
      end
      log_info "#{step_name} completed successfully".magenta.bold
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
