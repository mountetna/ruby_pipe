require 'fileutils'
module Pipeline
  module Script
    include Pipeline::Base
    include Pipeline::Logger
    include Pipeline::Usage

    module ClassMethods 
      attr_reader :steps
      def runs_steps(*step_list)
        @steps = step_list
      end
    end

    def self.included(base)
      base.extend(ClassMethods) if ClassMethods
    end

    class_var :steps

    def exclude_task? task
      nil
    end

    def initialize(d=nil)
      @defaults = d
    end

    def config
      @config ||= self.class.sister_class(:config).new(self,@defaults)
    end

    def run_action(*args)
      # run the current action
      send config.action, args
    end

    def schedule(args)
      # handle any previous steps
      curr_step = config.step

      create_step(config.prev_step).complete if config.prev_step

      return if config.single_step

      # schedule the execution for the current step
      job, splits = create_step(curr_step).setup_exec

      # schedule the scheduler for the next step
      create_step(config.next_step).setup_scheduler(job, splits) if config.next_step
    end

    def create_step(s)
      self.class.sister_class(s).new(self)
    end

    def exec(args)
      step = create_step(config.step)

      # clean up if you manage to get through the whole step.
      step.cleanup if step.exec 
    end

    def start_pipe(step)
      FileUtils.rm(config.error_file) if File.exists?(config.error_file)
      config.set_config :scheduler, scheduler_type
      job, splits = create_step(step).setup_exec
      create_step(config.next_step).setup_scheduler(job,splits) if config.next_step
    end

    usage "start <config_file.yml> [<step>]", "Start the pipeline at the beginning or at <step>"
    usage "run_step <config_file.yml> <step_name>", "Run just the named step"
    usage "audit <config_file.yml>", "Audit the pipeline to see which steps are complete."
    usage "stop [please]", "stop the pipeline, optionally waiting for the current step to finish"
    usage "generate <cohort_name>", "Generate a new config file for a cohort of samples."

    def generate(args)
      self.class.daughter_class(:config_generator).new *args
    end

    def run_step(args)
    end

    def start(args)
    end

    def init(args)
      # do something
      cmd = args.shift.to_sym

      if !usages[cmd]
        usage
        exit
      end

      ENV['CONFIG'] = args.shift if [ :start, :run_step, :audit, :stop ].include? cmd
      
      case cmd
      when :start
        config.set_config :action, :init
        start_pipe (args[0] || steps.first).to_sym
      when :run_step
        exit unless args.length == 1
        config.set_config :action, :init
        config.set_config :single_step, :true
        start_pipe args[0].to_sym
      when :audit
        config.set_config :action, :init
        audit(args)
      when :stop
        config.set_config :action, :init
        stop(args)
      when :generate
        generate(args)
      else
        usage
      end
    end
  end
end
