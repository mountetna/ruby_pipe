module Pipeline

  module Step
    def each_task(task)
      log_console "Files for #{config.splits} trials".yellow.bold

      files = {}
      if config.splits
        config.splits.times do |i|
          config.set_config :job_index, i
          task_files files
        end
      else
        config.set_config :job_index, 0
        task_files files
      end
      files.each do |t,f|
        next if task && t != task
        yield t, f
      end
    end

    def task_files mem
      self.class.tasks.each do |t|
        mem[t] ||= { :dump_files => [], :out_files => [] }
        task = create_task(t)
        mem[t][:dump_files].concat task.dump_files.map{|d| config.send d}.flatten
        mem[t][:out_files].concat task.out_files.map{|d| config.send d}.flatten
      end
    end
  end
  module Script
    def each_step step, task
      if step == :all
        steps.each do |s|
          create_step(s).each_task(task) do |t,f|
            yield t,f
          end
        end
      else
        abort "No such step" if !steps.include? step
        create_step(step).each_task(task) do |t,f|
          yield t,f
        end
      end
    end

    def file_size files
      files.map{|f| File.exists?(f) ? File.size(f) : 0}.reduce(:+) || 0
    end

    def file_summary files
      "#{files.count{|f| File.exists?f}} files, #{file_size(files).to_human}b"
    end

    def list step, task
      total_scratch = 0
      total_output = 0
      each_step step, task do |t,f|
        log_console "task #{t}".blue.bold
        log_console "  scratch files: #{file_summary f[:dump_files]}".red.bold
        log_console "  output files: #{file_summary f[:out_files]}".green.bold
        total_scratch += file_size f[:dump_files]
        total_output += file_size f[:out_files]
      end
      log_console "scratch total: #{total_scratch.to_human}b".cyan.bold
      log_console "output total: #{total_output.to_human}b".cyan.bold
    end

    def rm_files step, task, type
      total_size = 0
      total_count = 0
      each_step step, task do |t,f|
        f[type].each do |file|
          if File.exists? file
            total_size += File.size file
            total_count += 1
            if File.directory? file
              FileUtils.rm_rf file
            else
              File.unlink file 
            end
          end
        end
      end
      log_console "Deleted #{total_count} files (#{total_size.to_human}b)".cyan.bold
    end

    def scratch step, task
      rm_files step, task, :dump_files
    end

    def output step, task
      rm_files step, task, :out_files
    end

    def clean_job(args)
      cmd = args.shift.to_sym
      step = args.shift.to_sym

      task = args.first ? args.shift.to_sym : nil

      send cmd, step, task
    end
  end

  module Task
    module Cleaning
      class Cleaner
        def initialize(files,c)
          @files = files
          @config = c
        end
        def summary
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
      end

      def audit_files(files, c1, c2)
        aud = [ files.count do |f|
            filename = config.send(f)
            filename && File.size?(filename) && File.readable?(filename)
          end, files.size ]

        return aud + [ "(#{aud.first}/#{aud.last})".green ] #.bold.color( aud.first == aud.last ? c1 : c2 ) ]
      end

      def clean_list
        return if step.script.exclude_task? self
        log_console "  task #{task_name}:".blue.bold
        log_console "    dump files:"
        log_console "    out files:"
      end
    end
  end
end
