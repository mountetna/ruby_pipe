require 'fileutils'

module Pipeline
  module Step
    include Pipeline::Base
    include Pipeline::Logger
    include Pipeline::Scheduling
    # This creates a new step in a pipeline
    module ClassMethods
      attr_reader :run_chain, :tasks, :audit_vars
      def runs_task(*tasklist)
        @tasks = tasklist
      end
      alias :runs_tasks :runs_task

      def has_task *tasklist
        @available_tasks = tasklist
      end
      alias :has_tasks :has_task

      def available_tasks
        @available_tasks || tasks
      end

      def required
        @required ||= tasks.map{ |t| self.class.daughter_class(t).required_files }.flatten
      end

      def made
        @made ||= tasks.map{ |t| self.class.daughter_class(t).made_files }.flatten
      end

      def input
        required - made
      end

      def audit_report *vars
        @audit_vars = vars
      end

      def runs_on *job_item_keys
        @run_chain = job_item_keys
      end

      def resources(opts = nil)
        @resources ||= {:walltime => 3}
        @resources.update(opts || {})
      end
    end
    def self.included(base)
      base.extend(ClassMethods)
    end

    def init_hook
    end

    def initialize(s,tasks)
      @script = s

      init_hook

      set_tasks tasks if tasks.is_a? Array
      config.set_opt :step, step_name
      config.set_opt :job_array, job_array
      config.set_opt :resources, resources
      setup_logging unless config.action == :audit || config.action == :clean
    end

    class_var :run_chain, :tasks, :available_tasks, :resources, :audit_vars

    def set_tasks t
      @tasks = self.class.available_tasks.map{|e|
        t.include?(e) ? e : nil
      }.compact
    end

    def tasks
      @tasks ||= self.class.tasks
    end

    def job_array
      if run_chain
        obj = [ config ]
        run_chain.each do |key|
          obj = obj.collect(&key).flatten
        end
        return obj
      end
    end

    def step_name; self.class.name.split(/::/).last.snake_case.to_sym; end
    def script; @script; end
    def config; @script.config; end

    def setup_scheduler(prevjob, trials)
      # setup the scheduler to run this task.
      prevjob.strip!

      log_main "Scheduling #{step_name} with #{trials} trials".yellow.bold
      schedule_job :schedule, :wait => prevjob, :prev_trials => trials
    end

    def setup_exec
      # setup the scheduler to execute this task.

      log_main "Starting execution for #{step_name}".yellow.bold
      job = schedule_job :exec, { :trials => config.trials, :walltime => config.walltime }.merge(resources)
      [ job, config.trials ]
    end

    def make_error_file errors
      File.open(config.error_file, "w") do |f|
        f.puts "Script failed at:"
        errors.each do |step,error|
          f.puts "  #{step}:"
          f.puts error.map{|e| "    trial #{e[:trial]}:#{e[:task]}"}.join("\n")
        end
        f.puts "See logs in #{config.log_dir}/"
        f.puts "Resume with '#{config.pipe}_#{config.script} start #{config.config_file} #{config.step}'"
        f.puts "You may want to clean partial output of the last step before resuming, try '#{config.pipe}_#{config.script} clean #{config.config_file} list #{config.step}'"
        log_main "Script failed at #{config.step}".red.bold
      end
    end

    def get_errors
      errors = File.foreach(config.error_pid).map do |l|
        # if you were asked to stop, just exit, don't worry about errors
        exit if l.chomp == "stop"

        Hash[[ :step, :trial, :task ].zip( l.chomp.split)]
      end.group_by{ |l| l[:step] }
      FileUtils.rm(config.error_pid)
      return errors
    end

    def complete
      # consolidate logs for this step. 
      
      # Error out if there is an error pid
      if File.exists?(config.error_pid)
        make_error_file get_errors
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
      return # if config.keep_temp_files
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
      log_info "Starting #{step_name} trial #{config.job_index} on host #{config.hostname}".magenta.bold
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
