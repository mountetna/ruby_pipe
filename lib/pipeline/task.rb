module Pipeline
  module Task
    include Pipeline::Tools
    include Pipeline::Logger
    module ClassMethods
      @required_files = []
      @made_files = []
      def requires_file(*required)
        @required_files = required
      end
      alias :requires_files :requires_file
      def required_files; @required_files; end

      def makes_file(*made)
        @made_files = made
      end
      alias :makes_files :makes_file
      def made_files; @made_files; end
    end

    def self.included(base)
      base.extend(ClassMethods)
    end

    def initialize(s)
      @step = s
    end

    def config
      step.config
    end
    
    def step
      @step
    end

    def task_name
      self.class.name.split(/::/).last.snake_case
    end
    
    def made_files
      self.class.made_files
    end

    def required_files
      self.class.required_files
    end

    def should_run
      log_info "#{task_name}:"
      # this will exit the step if it is missing.
      check_required

      # this will merely move on to the next step
      return make_files?
    end

    def check_required
      required_files.each do |f|
        filename = config.send(f)
        error_exit "Could not find required file #{f} at #{filename}." if !File.size?(filename) 
        error_exit "File #{f} at #{filename} is not readable." if !File.readable?(filename)

        log_debug "#{f} ok at #{filename}" if config.verbose?
      end
    end

    def make_files?
      all_made = made_files.any?
      made_files.each do |f|
        filename = config.send(f)
        if File.size?(filename) && File.readable?(filename)
          log_debug "#{task_name}: #{f} already made at #{filename}" if config.verbose?
        else
          all_made = nil
        end
      end
      log_info "#{task_name}: skipping, all files made." if all_made
      return !all_made
    end

    def run
    end
  end
end
