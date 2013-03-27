require 'fileutils'

module Pipeline
  module Logger
    MAINLOG = File.open("/dev/null","w")

    def setup_logging
      FileUtils.mkdir_p config.log_dir
      STDOUT.reopen(config.step_log,"a")
      STDOUT.sync = true
      STDERR.reopen(STDOUT)
      STDERR.sync = true
      MAINLOG.reopen(config.main_log,"a")
    end

    def log_main(txt)
      log "PIPE", :blue, txt, MAINLOG
    end

    def log_info(txt)
      log "INFO", :blue, txt
    end

    def log_debug(txt)
      log "DEBUG", :red, txt.red
    end

    def log_error(txt)
      log "ERROR", :red, txt.red.bold
    end

    def log_output(txt)
      log "OUTPUT", :white, txt
    end

    def log_command(txt)
      log "COMMAND", :red, txt.cyan
    end

    def log_console(txt)
      s = ("%10s" % config.step).green.bold
      puts "%s %s" % [ s, txt ]
    end

    def log(t,color,txt,f=STDOUT)
      s = ("%10s" % config.step).green.bold
      t = ("%-10s" % t).send color
      date = ("%-20s" % Time.now.strftime("%Y-%m-%d %H:%M:%S")).white.bold
      f.puts "%s %s %s %s" % [ date, s, t, txt ]
    end
  end
end
