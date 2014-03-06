require 'pipeline'
require 'hash_table'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def_var :cohort_name do :blat end
  end
  class StartBlat
    include Pipeline::Step
    runs_task :start_blat

    def init_hook
      resources.update(:nodes => config.blat_node, :port => config.blat_port)
    end

    class StartBlat
      include Pipeline::Task
      no_requirements

      def run
        blat_start_server or error_exit "Couldn't start blat server"
      end
    end
  end

  class Blat
    include Pipeline::Script
    include Pipeline::Usage
    runs_steps :start_blat

    def start(args)
      args = set_config [File.join( ENV['LIB_DIR'], "config", "utility_blat.yml" )]
      start_pipe (args[0] || steps.first).to_sym
    end
    usage "start", "Start the blat server"

    def stop(args)
      args = set_config [File.join( ENV['LIB_DIR'], "config", "utility_blat.yml" )]
      config.set_config :scheduler, scheduler_type
      cancel_job
    end
    usage "stop", "Start the blat server"

    def_module :default, :start_blat => true
  end
end
