require 'time'
require 'pathname'
require 'hash_table'
require 'logger'

module Pipeline
  module Script
    include TimeFormat
    def read_log_file(file_loc)
      lines = IO.readlines(file_loc)
      filename = Pathname.new(file_loc).basename.to_s # get info from filename
      curr_step = filename.split('.')[2]
      curr_trial = filename.split('.')[3]
      if lines
        tasks = []
        curr_task = nil
        start_time = nil
        end_time = nil
        lines.each do |line|
          if /\x1b[^m]*m/.match(line) #only lines with ansi formatting should be checked
            l = line.gsub(/\x1b[^m]*m/, '') # remove ansi formatting
            date, time, step, type, text = l.split(" ",5)
            if text =~ /^task (\w+):/
              curr_task = $1
              start_time = "#{date} #{time}"
            end
            if text =~ /^#{curr_task} completed/
              end_time = "#{date} #{time}"
              tasks << {:step => curr_step,
                        :trial => curr_trial,
                        :task => curr_task,
                        :start_time => start_time,
                        :end_time => end_time,
                        :running_time => (Time.parse(end_time) - Time.parse(start_time)).to_s } #running time is in seconds
            end
          end
        end
      end
      tasks
    end

    def print_timer_metrics f
      o = HashTable.new nil
      o.header = @all_tasks.first.keys
      @all_tasks.each do |t|
        o.add_line t
      end
      o.print f
    end

    def analyze_logs(log_dir)
      logs = Dir.glob(log_dir + '/' + '*.log')
      @all_tasks = []
      logs.each do |log|
        tasks = read_log_file(log)
        if tasks.empty?
          puts "couldn't parse #{log}"
        else
          @all_tasks.concat tasks
        end
      end
    end

    private

    def mean_running_time task
      task.inject(0){|sum, t| sum += t[:running_time].to_i}/task.length.to_f
    end

    def show_task_stats task, step
      task = step.create_task task
      tasks = @all_tasks.select{|t| t[:step].to_sym == step.step_name && t[:task] == task.task_name}

      return if !tasks || tasks.empty?
      total_time = tasks.group_by{|t| t[:trial]}.map{|trial, tsk|
        mean_time = mean_running_time(tsk)
      }.inject(0){|sum,t| sum += t.to_i}
      mean_time = mean_running_time tasks
      max_time = tasks.map{|t| t[:running_time].to_i}.max
      puts "  #{task.task_name.blue.bold} mean: #{human_time(mean_time.to_i).red.bold} total: #{human_time(total_time).red.bold} max: #{human_time(max_time).red.bold}"
      [ mean_time, total_time, max_time ]
    end

    def show_step_stats step
      step = create_step step
      puts "#{step.step_name}".green.bold
      mean_time, total_time, max_time = [0,0,0]
      step.tasks.each do |task|
        stats = show_task_stats task, step
        next if !stats
        mean_time += stats[0]
        total_time += stats[1]
        max_time += stats[2]
      end
      puts "  mean: #{human_time(mean_time.to_i).magenta.bold} total: #{human_time(total_time).magenta.bold} max: #{human_time(max_time).magenta.bold}"
      [ mean_time, total_time, max_time ]
    end

    def show_script_stats
      mean_time, total_time, max_time = [0,0,0]
      steps.each do |step|
        stats = show_step_stats step
        next if !stats
        mean_time += stats[0]
        total_time += stats[1]
        max_time += stats[2]
      end
      puts "mean: #{human_time(mean_time.to_i).yellow.bold} total: #{human_time(total_time).yellow.bold} max: #{human_time(max_time).yellow.bold}"
    end

    def timer_job(args)
      analyze_logs config.log_dir
      #print_timer_metrics config.timer_metrics

      show_script_stats
    end
  end
end
