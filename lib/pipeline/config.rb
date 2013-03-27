require 'yaml'
module Pipeline
  module Config
    def def_var name, &block
      class_procs.update name => block
    end
    def class_procs; @class_procs ||= {}; end
  end
  module BaseConfig
    extend Pipeline::Config

    def_var :next_step do @script.steps[ @script.steps.index(step) + 1 ] end
    def_var :prev_step do @script.steps[ @script.steps.index(step) - @script.steps.size - 1 ] end
    def_var :pipe do @script.class.name.sub(/::.*/,"").snake_case.to_sym end
    def_var :script do @script.class.name.snake_case.to_sym end
    def_var :pipe_script do "#{pipe}_#{script}" end
    def_var :step_log do "#{log_dir}/#{pipe}.#{cohort_name}.#{step}.#{job_index}.log" end
    def_var :main_log do "#{log_dir}/#{pipe}.#{cohort_name}.#{script}.log" end

    def_var :error_pid do "#{cohort_scratch}/error.pid" end
    def_var :error_file do "ERROR.#{cohort_name}" end
    def_var :verbose? do verbose == "yes" || verbose == "true" end
    def_var :config_dir do "#{lib_dir}/config" end
    def_var :tools_config do "#{config_dir}/tools.yml" end

    def_var :job_index do job_number ? job_number - 1 : 0 end
    def_var :job_item do job_array[job_index] end
    def_var :threads do resources[:threads] end
    def_var :qual_type do "phred64" end

    def_var :cohort_dir do |dir,s| File.join dir, cohort_name end
    def_var :cohort_scratch do cohort_dir scratch_dir end
    def_var :cohort_scratch_file do |affix| File.join cohort_scratch, affix end 
    def_var :cohort_output do cohort_dir output_dir end
    def_var :cohort_output_file do |affix| File.join cohort_output, affix end 
    def_var :cohort_metrics do cohort_dir metrics_dir end
    def_var :cohort_metrics_file do |affix| File.join metrics_dir, "#{cohort_name}.#{affix}" end 

    def_var :sample_dir do |dir,s| File.join dir, (s || sample_name) end
    def_var :sample_scratch do |s| sample_dir scratch_dir, s end
    def_var :sample_scratch_file do |affix,s| File.join sample_scratch(s), affix end 
    def_var :sample_output do |s| sample_dir output_dir, s end
    def_var :sample_output_file do |affix,s| File.join sample_output(s), affix end 
    def_var :sample_metrics do |s| sample_dir metrics_dir, s end
    def_var :sample_metrics_base do |s| File.join metrics_dir, "#{s || sample_name}" end 
    def_var :sample_metrics_file do |affix,s| File.join metrics_dir, "#{s || sample_name}.#{affix}" end 
    def_var :sample_name do sample.sample_name end

    def_var :reference_fa do hg19_fa end
    def_var :reference_name do :hg19 end
    def_var :reference_date do :feb_2009 end
    def_var :reference_gtf do hg19_ucsc_gtf end


    def splits
      job_array ? job_array.size : nil
    end

    def sample(name=nil)
      return samples.find { |s| s[:sample_name] == name } if name
      job_item.parent_with_property :sample_name if job_array
    end

    def tumor_samples
      n = samples.select { |s| s[:normal_name] }
      n.empty? ? samples[1..-1] : n
    end

    def chroms
      # read chromosomes from the fasta dict
      @chroms ||= File.foreach("/taylorlab/resources/human/hg19/ucsc_feb_2009/hg19.dict").map do |s|
        s.match(/@SQ.*SN:(\w+)\s/) do |m| 
          next if m[1] == "chrM"
          { :contig => m[1] }
        end
      end.compact
    end

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

      @config = {}

      # also include the script config
      #self.class.send(:include, script.class.daughter_class(:config))

      load_procs
      load_env_vars
      load_tools_config

      set_dir

      load_config
    end

    def load_procs
      @procs = {}
      self.class.ancestors.each do |a|
        @procs.update a.class_procs if a.respond_to? :class_procs
      end
    end

    def load_env_vars
      # if these exist, put them in @opts
      Hash[
        :PBS_O_WORKDIR => :work_dir,
        :ACTION => :action,
        :STEP => [ :step, :to_sym ],
        :LIB_DIR => :lib_dir,
        :MOAB_JOBARRAYINDEX => [:job_number, :to_i],
        :PBS_ARRAYID => [:job_number, :to_i],
        :MOAB_JOBARRAYRANGE => [:array_range, :to_i ],
        :SINGLE_STEP => :single_step,
        :SCHEDULER => :scheduler,
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
      @opts.update( YAML.load_file( tools_config ) )

      # update the path
      ENV['PATH'] = ( tools_path.map{ |t| send t.to_sym } + ENV['PATH'].split(/:/) ).join ":" 
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
      @config=Pipeline::SampleObject.new YAML.load_file(config_file), self
    end

    def method_missing(meth,*args,&block)
      # always look in the config file first, so we can overrule things there
      meth = meth.to_sym if meth.is_a? String
      if @config[meth]
        return @config[meth]
      # if not there, it could be in the procs table, to be found interactively
      elsif @procs[meth]
        return instance_exec( *args, &@procs[meth])
      # if not there, it is likely in options. These can be nil for empty variables
      elsif @opts.has_key? meth
        return @opts[meth]
      else
        super
      end
    end
  end
end
