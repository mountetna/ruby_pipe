require 'colored'

module Pipeline
  module Logger
    MAINLOG = File.open("/dev/null","w")

    def setup_logging
      STDOUT.reopen(config.step_log,"a")
      STDOUT.sync = true
      STDERR.reopen(STDOUT)
      STDERR.sync = true
      MAINLOG.reopen(config.main_log,"a")
    end

    def log_info(txt)
      log "INFO", :blue, config.step, txt
    end

    def log_debug(txt)
      log "DEBUG", :red, config.step, txt.red
    end

    def log_error(txt)
      log "ERROR", :red, config.step, txt.red.bold
    end

    def log_output(txt)
      log "OUTPUT", :white, config.step, txt
    end

    def log(t,color,s,txt)
      s = ("%10s" % s).green.bold
      t = ("%-10s" % t).send color
      date = ("%-20s" % Time.now.strftime("%Y-%m-%d %H:%M:%S")).white.bold
      puts "%s %s %s %s" % [ date, s, t, txt ]
    end
  end
end
