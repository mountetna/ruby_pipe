module Pipeline
  class Task
    class << self
      # list of required filenames, strings, etc.
      # that we need as input

      def input name, options={}
        # we need to know what name is to run
        require_symbol(name, options.merge(type: :input))
      end

      def output name, options={}
        require_symbol(name, options.merge(type: :output))
      end

      def resource name
        @resources ||= []
        @resources.push name
      end

      def tool name
        @tools ||= []
        @tools.push name
      end

      private

      def require_symbol name, options)
        @required ||= {}
        @required[name] = options)
      end
    end

    def initialize(config)
      @inputs = []
      @outputs = []
      @resources = []
      @tools = []

      create_requirements(config)
    end

    def create_requirements config
      requirements = self.class.required

      missing = requirements.reject do |name, opts|
        config.has_key?(name) || opts[:optional]
      end

      raise Pipeline::Error, "Missing #{missing.keys} as input to #{task_name}!" unless missing.empty?
    end

    def should_run?
      # You are not marked as 'done' or 'running'
      !(running? || done?) &&

      # You have not made the correct output
      !satisfied? &&

      # You have the correct inputs
      input_exists? &&

      # You have the correct tools
      tools_exist? &&

      # You have the correct resources
      resources_exist?
    end

    private

    def error_exit(txt)
      raise Pipeline::PipeError, txt
    end

    def error_out(txt)
      exit
    end

    def has_required?
    end

    def make_files?
      all_made = made_files.any?
      made_files.each do |f,type|
        filename = config.send(f)
        if filename.is_a? Array
          filename.each { |fn| all_made = nil if make_file? fn, f, type }
        else
          all_made = nil if make_file? filename, f, type
        end
      end
      log_info "#{task_name}: skipping, all files made".blue if all_made
      return !all_made
    end

    def error_check
      begin
        yield
      rescue Pipeline::PipeError => e
        error_out e.message
      rescue => e
        # there was a ruby error, bail out.
        puts e, e.backtrace
        error_out "script error"
      end
    end

    def exec
      error_check do
        run
      end
      log_info "#{task_name} completed successfully".green.bold
    end

    def run
    end
  end
end
