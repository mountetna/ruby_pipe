require 'yaml'
module Pipeline
  module Config
    # This creates a new step in a pipeline
    def initialize(script,opts=nil)
      @script = script
      # stuff that should default to nil goes here, rather than using method_missing
      @opts = {
        :work_dir => nil,
        :verbose => nil,
        :job_number => nil,
        :keep_temp_files => nil,
        :single_step => nil
      }.merge(opts || {})
      @procs = {}
      @config = {}

      load_env_vars

      load_tools_config

      # also include the script config
      self.class.send(:include, script.class.daughter_class(:config))

      set_dir

      load_config

      load_procs
    end

    def load_env_vars
      # if these exist, put them in @opts
      Hash[
        :PBS_O_WORKDIR => :work_dir,
        :ACTION => :action,
        :STEP => [ :step, :to_sym ],
        :LIB_DIR => :lib_dir,
        :MOAB_JOBARRAYINDEX => [:job_number, :to_i],
        :PBS_ARRAYINDEX => [:job_number, :to_i],
        :MOAB_JOBARRAYRANGE => [:array_range, :to_i ],
        :SINGLE_STEP => :single_step,
        :CONFIG => :config_file
      ].each do |ev,o|
        if o.is_a? Array
          conv = o.last; o = o.first
        end
        next if !ENV[ev.to_s] || @opts[o]
        if conv
          @opts[o] = ENV[ev.to_s].send conv
        else
          @opts[o] = ENV[ev.to_s]
        end
      end
    end

    def load_tools_config
      @opts.update( YAML.load_file( "/taylorlab/resources/tools.yml" ) )
    end

    def set_opt(o,v)
      @opts[o] = v
    end

    def set_config(c,v)
      @config[c] = v
    end

    def set_dir
      Dir.chdir(work_dir) if work_dir
    end

    def load_config
      @config=YAML.load_file(config_file)
    end

    def load_procs
      @procs = {
        :next_step => proc { @script.steps[ @script.steps.index(step) + 1 ] },
        :prev_step => proc { @script.steps[ @script.steps.index(step) - @script.steps.size - 1 ] },
        :pipe => proc { self.class.name.sub(/::.*/,"").to_sym },
        :script => proc { @script.class.name.snake_case.to_sym },
        :step_log => proc { "#{log_dir}/#{pipe}.#{job_name}.#{step}.#{job_index}.log" },
        :main_log => proc { "#{log_dir}/#{pipe}.#{job_name}.#{script}.log" },
        :job_scratch => proc { "#{scratch_dir}/#{job_name}" },
        :sample_scratch => proc { "#{scratch_dir}/#{sample_name}" },
        :error_pid => proc { "#{job_scratch}/error.pid" },
        :error_file => proc { "ERROR.#{job_name}" },
        :verbose? => proc { verbose == "yes" || verbose == "true" },
        :job_index => proc { job_number ? job_number - 1 : 0 }
      }
    end

    def method_missing(meth,*args,&block)
      # always look in the config file first, so we can overrule things there
      meth = meth.to_sym if meth.is_a? String
      if @config[meth]
        return @config[meth]
      # if not there, it could be in the procs table, to be found interactively
      elsif @procs[meth]
        return @procs[meth].call(*args)
      # if not there, it is likely in options. These can be nil for empty variables
      elsif @opts.has_key? meth
        return @opts[meth]
      else
        super
      end
    end
  end
end
