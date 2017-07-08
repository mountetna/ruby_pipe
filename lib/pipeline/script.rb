require 'fileutils'
module Pipeline
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

      def across *lists, &block
      end

      def task name, reqs
        @tasks ||= {}
        @tasks[name] = reqs
      end
    end
  end
end
