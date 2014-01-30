require 'fileutils'
module Pipeline
  module Task
    module ClassMethods
      def class_init
        @required_files = []
        @skip_symbols = []
        @dump_files = []
        @out_files = []
      end

      def skip_without *required
        requires_file *required
        @skip_symbols = required
      end

      def requires_file(*required)
        @required_files = required
      end
      alias :requires_files :requires_file

      def no_requirements
        @required_files = []
      end

      def dumps_file(*dump)
        @dump_files = dump
      end
      alias :dumps_files :dumps_file

      def outs_file(*outs)
        @out_files = outs
      end
      alias :outs_files :outs_file

      def made_files; @out_files + @dump_files; end
      attr_reader :required_files, :skip_symbols, :dump_files, :out_files
    end
  end
  class PipeError < StandardError
  end
  module Task
    module ErrorHandling
      def error_exit(txt)
        raise Pipeline::PipeError, txt
      end

      def error_out(txt)
        log_error txt

        log_error "Exiting at #{task_name}"

        log_main "#{config.step} trial #{config.job_index} failed to complete at #{task_name}.".red.bold

        open_error_pid do |f|
          f.puts "#{config.step} #{config.job_index} #{task_name}"
        end
        exit
      end
    end
  end

  module Task
    include Pipeline::Tools
    include Pipeline::Logger
    include Pipeline::Task::ErrorHandling
    include Pipeline::Task::Auditing
    include Pipeline::Task::Cleaning
    include Pipeline::Base

    def self.included(base)
      base.extend(ClassMethods)
      base.class_init if base.respond_to? :class_init
    end

    def initialize(s)
      @step = s
    end

    def config
      step.config
    end
    
    attr_reader :step

    def task_name
      self.class.name.split(/::/).last.snake_case
    end
    
    class_var :dump_files, :out_files, :made_files, :required_files, :skip_symbols

    def should_run
      # don't even log it if the script says to skip it
      error_check do
        return if step.script.exclude_task? self

        log_info "task #{task_name}:".white.bold

        # if all of the files exist, you shouldn't run
        return nil if !make_files?

        # this will exit the step if it is missing files.
        # it will skip if some required symbols are nil
        return nil if !check_required
        # this will merely move on to the next step
      end
      true
    end

    def check_file(filename,f)
      error_exit "Could not get filename for #{f}" if !filename

      # make sure the directory is there
      ensure_dir filename

      error_exit "Could not find required file #{f} at #{filename}" if !File.size?(filename) 
      error_exit "File #{f} at #{filename} is not readable" if !File.readable?(filename)

      log_debug "#{f} ok at #{filename}" if config.verbose?
    end

    def check_required
      skip_symbols.each do |s|
        return nil if !config.send s
      end

      required_files.each do |f|
        filename = config.send(f)
        if filename.is_a? Array
          filename.each {|ff| check_file ff, f }
        else
          check_file filename, f
        end
      end
      true
    end

    def make_file?(filename,f)
      if File.size?(filename) && File.readable?(filename)
        log_debug "#{task_name}: #{f} already made at #{filename}" if config.verbose?
        nil
      else
        ensure_dir filename
        true
      end
    end

    def make_files?
      all_made = made_files.any?
      made_files.each do |f|
        filename = config.send(f)
        if filename.is_a? Array
          filename.each { |fn| all_made = nil if make_file? fn, f }
        else
          all_made = nil if make_file? filename, f
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
