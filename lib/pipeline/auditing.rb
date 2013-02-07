module Pipeline
  module Script
    def audit(args)
      config.set_config :action, :audit
      steps.each do |s|
        step = create_step s

        config.splits.times do |i|
          config.set_config :job_index, i
          step.audit
        end
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
end
