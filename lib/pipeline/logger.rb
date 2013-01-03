require 'colored'

module Pipeline
  module Logger
    MAINLOG = File.open("/dev/null","w")

    def setup_logging
      STDOUT.reopen(config.step_log,"a")
      STDERR.reopen(STDOUT)
      MAINLOG.reopen(config.main_log,"a")
    end

    def log_info(txt)
      log "INFO", config.step, txt
    end

    def log_debug(txt)
      log "DEBUG", config.step, txt.red
    end

    def log_error(txt)
      log "ERROR", config.step, txt.red.bold
    end

    def log_output(txt)
      log "OUTPUT", config.step, txt
    end

    def log(t,s,txt)
      s = ("%10s" % s).green.bold
      t = ("%-10s" % t).red
      date = ("%-20s" % Time.now.strftime("%Y-%m-%d %H:%M:%S")).white.bold
      puts "%s %s %s %s" % [ date, s, t, txt ]
    end
  end
end
