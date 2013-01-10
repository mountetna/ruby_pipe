module Pipeline
  module Script
    include Pipeline::Logger
    module ClassMethods 
      def runs_steps(*step_list)
        @steps = step_list
      end

      def steps
        @steps
      end
    end

    def self.included(base)
      base.extend(ClassMethods)
    end

    def steps
      self.class.steps
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
      create_step(config.prev_step).complete if config.prev_step

      return if config.single_step

      # schedule the execution for the current step
      job = create_step(config.step).setup_exec

      # schedule the scheduler for the next step
      create_step(config.next_step).setup_scheduler(job) if config.next_step
    end

    def exec(args)
      step = create_step(config.step)

      # clean up if you manage to get through the whole step.
      step.cleanup if step.exec 
    end

    def start_pipe(step)
      FileUtils.rm(config.error_file) if File.exists?(config.error_file)
      job = create_step(step).setup_exec
      create_step(config.next_step).setup_scheduler(job) if config.next_step
    end

    def init(args)
      # do something
      cmd = args.shift

      ENV['CONFIG'] = args.shift if [ "start", "run_step", "audit" ].include? cmd
      
      case cmd
      when "start"
        start_pipe (args[0] || steps.first).to_sym
      when "run_step"
        exit unless args.length == 1
        config.set_config :single_step, :true
        start_pipe args[0].to_sym
      when "audit"
        audit(args)
      end
    end

    def audit(args)
      steps.each do |s|
        step = create_step s
        # config thinks we're now on this step
        #config.set_config :array_range, config.splits
        config.splits.times do |i|
          config.set_config :job_index, i
          step.audit
        end
      end
    end

    def create_step(s)
      self.class.sister_class(s).new(self)
    end
  end
end
