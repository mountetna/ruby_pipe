module Pipeline
  class Script
    def self.steps(*step_list)
      @@steps = step_list
    end

    def initialize
    end

    def config
      @config ||= self.class.sister_class(:config).new(self)
    end

    def run_action
      # run the current action
      send config.action
    end

    def schedule_action
      # handle any previous steps
      create_step(config.prev_step).cleanup if config.prev_step

      # schedule the execution for the current step
      job = create_step(config.step).setup_exec

      # schedule the scheduler for the next step
      create_step(config.next_step).setup_scheduler(job) if config.next_step
    end

    def exec_action
      step = create_step(config.step)

      # clean up if you manage to get through the whole step.
      step.cleanup if step.exec 
    end

    def create_step(s)
      self.class.sister_class(s).new(self)
    end
  end
end
