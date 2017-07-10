module Pipeline
  class TaskBlock
    def initialize list_names, &block
      @list_names = list_names
      @block = block
    end
  end

  class Script
    class << self
      def scratch name, file_loc
        @scratch ||= {}
        @scratch[name] = file_loc
      end

      def output name, file_loc
        @output ||= {}
        @output[name] = file_loc
      end

      def across *list_names, &block
        @task_blocks ||= []
        @task_blocks.push TaskBlock.new(list_names,&block)
      end
    end

    def tasks
      @tasks ||= []
      self.class.task_blocks.each do |task_block|
        lists = get_lists(task_block.list_names)
        instance_exec(*lists,&task_block.block)
      end
    end

    def task name, reqs
      @tasks ||= {}
      @tasks[name] = get_task(reqs)
    end
  end
end
