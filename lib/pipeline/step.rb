require 'fileutils'

module Pipeline
  module Step
    include Pipeline::Base
    include Pipeline::Logger
    include Pipeline::Scheduling
    # This creates a new step in a pipeline
    module ClassMethods
      attr_reader :job_items, :tasks, :available_tasks
      def runs_task(*tasklist)
        @tasks = tasklist
      end
      alias :runs_tasks :runs_task

      def has_task *tasklist
        @available_tasks = tasklist
      end
      alias :has_tasks :has_task

      def required
        @required ||= @tasks.map{ |t| self.class.daughter_class(t).required_files }.flatten
      end
      def made
        @made ||= @tasks.map{ |t| self.class.daughter_class(t).made_files }.flatten
      end
      def input
        required - made
      end

      def job_list &block
        @job_array = block
      end

      def runs_on *job_item_keys
        @job_items = job_item_keys
      end

      def resources(opts = nil)
        @resources ||= {:walltime => 3}
        @resources.update(opts || {})
      end
    end
    def self.included(base)
      base.extend(ClassMethods)
    end

    def initialize(s)
      @script = s
      config.set_opt :step, step_name
      config.set_opt :job_array, job_array
      config.set_opt :resources, resources
      setup_logging unless config.action == :audit || config.action == :clean
    end

    class_var :job_items, :tasks, :available_tasks, :resources

    def job_array
      if job_items
        obj = config.send job_items.first
        if job_items.size > 1
          obj = obj.collect(&job_items[1]).flatten
        end
        return obj
      end
    end

    def step_name; self.class.name.split(/::/).last.snake_case.to_sym; end
    def script; @script; end
    def config; @script.config; end

    def setup_scheduler(prevjob, splits)
      # setup the scheduler to run this task.
      prevjob.strip!

      log_main "Scheduling #{step_name}".yellow.bold
      schedule_job :schedule, :wait => prevjob, :prev_splits => splits
    end

    def setup_exec
      # setup the scheduler to execute this task.

      log_main "Starting execution for #{step_name}".yellow.bold
      job = schedule_job :exec, :splits => config.splits, :walltime => resources[:walltime], :threads => resources[:threads]
      [ job, config.splits ]
    end

    def complete
      # consolidate logs for this step. 
      
      # Error out if there is an error pid
      if File.exists?(config.error_pid)
        FileUtils.rm(config.error_pid)
        File.open(config.error_file, "w") do |f|
          f.puts "Script failed at #{config.step}, see logs in #{config.log_dir}/. Resume with '#{config.pipe}_#{config.script} start #{config.config_file} #{config.step}'"
          log_main "Script failed at #{config.step}".red.bold
        end
        exit
      end

      # You're good, let them know.
      log_main "Step #{config.step} completed successfully.".green.bold
    end

    def vacuum
    end

    def cleanup
      # run the vacuum script
      vacuum
      return if config.keep_temp_files
      tasks.each do |t|
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
      log_info "Starting #{step_name} trial #{config.job_index}".magenta.bold
      tasks.each do |t|
        task = create_task t

        task.exec if task.should_run
      end
      log_info "#{step_name} completed trial #{config.job_index} successfully".magenta.bold
      log_main "#{step_name} completed trial #{config.job_index} successfully".magenta.bold
      return true
    end

    def create_task(t)
      self.class.daughter_class(t.to_sym).new(self)
    end
  end
end
