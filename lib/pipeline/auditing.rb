module Pipeline
  module Script
    def audit_step s
      step = create_step s
      log_console "runs on #{(step.job_items || [:cohort]).join(".")}".cyan.bold

      if config.splits
        config.splits.times do |i|
          config.set_opt :job_number, i+1
          step.audit
        end
      else
        config.set_opt :job_number, nil
        step.audit
      end
    end

    def audit_job(args)
      if args.last == "verbose"
        args.pop
        config.set_config :verbose, true
      end
        
      if args.first
        step = args.first.to_sym
        abort "No such step!".red.bold if !steps.include? step
        audit_step step
      else
        steps.each do |s|
          audit_step s
        end
      end
    end
  end

  module Step
    def audit
      log_console "Auditing #{self.class.name.snake_case} trial #{config.job_index}".yellow.bold

      if self.audit_vars
        self.audit_vars.each do |v|
          val = config.send(v) || "nil".red
          log_console "#{v}: ".magenta.bold + val
        end
      end

      tasks.each do |t|
        create_task(t).audit
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

        def symbols_missing?
          symbol_count != total
        end

        def symbol_count
          @files.map{|f| @config.send f}.compact.size
        end

        def needed?
          total != 0 && missing?
        end

        def count
          @count ||= @files.count do |f|
            begin
              filename = @config.send f
            rescue => e
              puts e, e.backtrace
              puts "Could not successfully retrieve symbol :#{f}".red.bold
              exit
            end
            if filename && filename.is_a?(Array)
              filename.all? do |fn|
                fn && File.size?(fn) && File.readable?(fn)
              end
            else
              filename && File.size?(filename) && File.readable?(filename)
            end
          end
        end
        def print_files no_array=nil
          @files.each do |f|
            filename = @config.send f
            if filename && filename.is_a?(Array) && !no_array
              filename.each do |fn|
                yield(f,fn)
              end
            else
              yield(f,filename)
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

      def empty? fn
        !fn || !File.size?(fn) || !File.readable?(fn)
      end

      def audit
        return if step.script.exclude_task? self
        log_console "task #{task_name}:".blue.bold
        sym = Audit.new(skip_symbols,config)
        req = Audit.new(required_files,config)
        dump = Audit.new(dump_files,config)
        out = Audit.new(out_files,config)
        log_console "required files: #{req.summary}" if req.total > 0
        log_console "dump files: #{dump.summary}" if dump.total > 0
        log_console "out files: #{out.summary}" if out.total > 0
        if sym.symbols_missing?
          log_console "doesn't need to run".red
          sym.print_files do |f,fn| 
            log_console "#{f} => #{fn}".blue if !fn
          end
          return
        end
        if !dump.missing? && !out.missing?
          log_console "won't run, all files are present".blue
          out.print_files do |f,fn|
            log_console "#{f} => #{fn} (#{File.human_size(fn)})".green
          end if config.verbose
          return
        end
        if req.missing?
          log_console "won't run, required files are missing".red
          req.print_files do |f,fn| 
            log_console "#{f} => #{fn}".red if empty? fn
          end
          return
        end
        if dump.missing?
          log_console "will run, needs to make dump files".green
          dump.print_files do |f,fn| 
            log_console "#{f} => #{fn}".red if empty? fn
          end
        end
        if out.missing?
          log_console "will run, needs to make out files".green
          out.print_files do |f,fn|
            log_console "#{f} => #{fn}".red if empty? fn
          end
        end
      end
    end
  end
end
