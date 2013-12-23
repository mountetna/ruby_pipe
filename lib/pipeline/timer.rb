require 'time'
require 'pathname'
require 'hash_table'
require 'logger'

module Pipeline
  module Script
    TASKS = []

    def read_log_file(file_loc)
      lines = IO.readlines(file_loc)
      filename = Pathname.new(file_loc).basename.to_s # get info from filename
      curr_step = filename.split('.')[2]
      curr_trial = filename.split('.')[3]
      if not lines.nil?
        tasks = []
        curr_task = nil
        start_time = nil
        end_time = nil
        lines.each do |line|
          if /\x1b[^m]*m/.match(line) #only lines with ansi formatting should be checked
            l = line.gsub(/\x1b[^m]*m/, '') # remove ansi formatting
            l = l.split(" ") #ex: ["2013-11-19", "05:12:30", "align", "INFO", "task", "mark_duplicates:"]
            if l[4] == "task"
              curr_task = l[5].gsub(":", "")
              start_time = "#{l[0]} #{l[1]}"
            end
            if l[4] == curr_task and l[5].include? "completed"
              end_time = "#{l[0]} #{l[1]}"
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

    def output (f, tasks)
      o = HashTable.new ''
      o.header = tasks[0].keys
      tasks.each do |t|
        o.add_line t
      end
      o.print f
    end

    def loop_through_files(log_dir)
      logs = Dir.glob(log_dir + '/' + '*.log')
      logs.each do |log|
        tasks = read_log_file(log)
        if not tasks.empty?
          TASKS << tasks
        else
          puts "couldn't parse #{log}"
        end
      end
    end

    def main()
      loop_through_files("./log")
      output("./metrics/timer_metrics", TASKS.flatten)
    end

    def timer_job(args)
      main()
    end
  end
end

