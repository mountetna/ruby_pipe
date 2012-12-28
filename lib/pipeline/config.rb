module Pipeline
  module Config
    # This creates a new step in a pipeline
    def initialize(script)
      @script = script
      # read in the config file
      set_pbs_dir

      load_config

      CONFIG={}
    end

    def set_pbs_dir
      Dir.chdir(ENV['PBS_O_WORKDIR']) if ENV['PBS_O_WORKDIR']
    end

    def load_config
      File.foreach(ENV['CONFIG']) do |f|
        next if f =~ /^#/
        f.match( /\s*([\w]+)=(.*)$/ ) do |m|
          CONFIG[m[1]] = m[2]
        end
      end
    end

    def action
      (ENV['ACTION'] + "_action").to_sym
    end

    def step
      ENV['STEP'].to_sym
    end
  end
end
