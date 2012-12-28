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
    end
    def self.included(base)
      base.extend(ClassMethods)
    end

    def initialize(s)
      @script = s
    end

    def script
      @script
    end

    def config
      @script.config
    end

    def step
      self
    end

    def schedule_schedule(prevjob)
      # setup the scheduler to run this task.
    end

    def exec
      # execute this step. Find tasks that have not yet been completed.
      self.class.tasks.each do |t|
        task = create_task t

        # execute the task
        task.run if task.should_run
      end
      return true
    end

    def step_name
      self.class.name.split(/::/).last
    end

    def create_task(t)
      self.class.daughter_class(t).new(self)
    end
  end
end
