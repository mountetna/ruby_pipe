require 'fileutils'
module Pipeline
  module Task
    module ClassMethods
      def class_init
        @required_files = []
        @dump_files = []
        @out_files = []
      end

      def requires_file(*required)
        @required_files = required
      end
      alias :requires_files :requires_file
      def required_files; @required_files; end

      def dumps_file(*dump)
        @dump_files = dump
      end
      alias :dumps_files :dumps_file
      def dump_files; @dump_files; end

      def outs_file(*outs)
        @out_files = outs
      end
      alias :outs_files :outs_file
      def out_files; @out_files; end

      def made_files; @out_files + @dump_files; end
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

        ensure_dir config.job_scratch

        File.open(config.error_pid,"w") do |f|
          f.puts "failed"
        end
        exit
      end
    end
  end

  module Task
    module Auditing
      class Audit
        def initialize(files,c)
          @files = files
          @config = c
        end
        def summary
          "(#{count}/#{total})".send( count == total ? :green : :red )
        end

        def total
          @files.size
        end

        def missing?
          count != total
        end

        def needed?
          total != 0 && missing?
        end

        def count
          @count ||= @files.count do |f|
            filename = @config.send f
            if filename && filename.is_a?(Array)
              filename.all? do |fn|
                fn && File.size?(fn) && File.readable?(fn)
              end
            else
              filename && File.size?(filename) && File.readable?(filename)
            end
          end
        end
        def print_missing
          @files.each do |f|
            filename = @config.send f
            if filename && filename.is_a?(Array)
              filename.each do |fn|
                yield(f,fn) if !fn || !File.size?(fn) || !File.readable?(fn)
              end
            else
              yield(f,filename) if !filename || !File.size?(filename) || !File.readable?(filename)
            end
          end
        end
      end

      def audit_files(files, c1, c2)
        aud = [ files.count do |f|
            filename = config.send(f)
            filename && File.size?(filename) && File.readable?(filename)
          end, files.size ]

        return aud + [ "(#{aud.first}/#{aud.last})".green ] #.bold.color( aud.first == aud.last ? c1 : c2 ) ]
      end

      def audit
        return if step.script.exclude_task? self
        log_info "task #{task_name}:".blue.bold
        req = Audit.new(required_files,config)
        dump = Audit.new(dump_files,config)
        out = Audit.new(out_files,config)
        log_info "required files: #{req.summary}" if req.total > 0
        log_info "dump files: #{dump.summary}" if dump.total > 0
        log_info "out files: #{out.summary}" if out.total > 0
        if req.missing?
          log_info "won't run, required files are missing".red
          req.print_missing { |f, filename| log_info "#{f} => #{filename}".red }
          return
        end
        if dump.missing?
          log_info "will run, needs to make dump files".green
          dump.print_missing { |f, filename| log_info "#{f} => #{filename}".red }
        end
        if out.missing?
          log_info "will run, needs to make out files".green
          out.print_missing { |f, filename| log_info "#{f} => #{filename}".red }
        end
        if !dump.missing? && !out.missing?
          log_info "won't run, all files are present".blue
        end
      end
    end
  end

  module Task
    include Pipeline::Tools
    include Pipeline::Logger
    include Pipeline::Task::ErrorHandling
    include Pipeline::Task::Auditing

    def self.included(base)
      base.extend(ClassMethods)
      base.class_init if base.respond_to? :class_init
    end

    def initialize(s)
      @step = s
    end

    def ensure_dir(*dirs)
      dirs.each do |f|
        FileUtils.mkdir_p f
      end
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
    
    def dump_files; self.class.dump_files; end
    def out_files; self.class.out_files; end
    def made_files; self.class.made_files; end
    def required_files; self.class.required_files; end

    def should_run
      # don't even log it if the script says to skip it
      return if step.script.exclude_task? self

      log_info "task #{task_name}:".white.bold

      # this will exit the step if it is missing.
      check_required

      # this will merely move on to the next step
      return make_files?
    end

    def check_file(filename,f)
      error_exit "Could not get filename for #{f}" if !filename

      # make sure the directory is there
      ensure_dir File.dirname(filename)

      error_exit "Could not find required file #{f} at #{filename}" if !File.size?(filename) 
      error_exit "File #{f} at #{filename} is not readable" if !File.readable?(filename)

      log_debug "#{f} ok at #{filename}" if config.verbose?
    end

    def check_required
      required_files.each do |f|
        filename = config.send(f)
        if filename.is_a? Array
          filename.each {|ff| check_file ff, f }
        else
          check_file filename, f
        end
      end
    end

    def make_file?(filename,f)
      if File.size?(filename) && File.readable?(filename)
        log_debug "#{task_name}: #{f} already made at #{filename}" if config.verbose?
        nil
      else
        ensure_dir File.dirname(filename)
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

    def exec
      begin
        run
      rescue Pipeline::PipeError => e
        error_out e
      rescue => e
        # there was a ruby error, bail out.
        puts e, e.backtrace
        error_out "script error"
      end
      log_info "#{task_name} completed successfully".green.bold
    end

    def run
    end
  end
end
