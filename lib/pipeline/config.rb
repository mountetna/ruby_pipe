require 'yaml'
module Pipeline
  module Config
    def def_var name, &block
      class_procs.update name => block
    end
    def empty_var *symbols
      class_opts.update(symbols.map{|s| { s => nil } }.reduce :merge)
    end
    def class_opts; @class_opts ||= {}; end
    def class_procs; @class_procs ||= {}; end
    def dir_tree tree
      @the_tree ||= {}
      @the_tree.deep_merge! tree
      get_procs_from_tree tree
    end

    def job_items *jis
      jis.each do |ji|
        self.instance_eval <<-EOT
          def_var :#{ji} do |item| find_job_item :#{ji}, item; end
        EOT
      end
    end

    def get_procs_from_tree tree, ancestor=nil
      tree.each do |file,blob|
        full_name = [ ancestor, file ].compact.join("/")
        case blob
        when Hash
          get_procs_from_tree blob, full_name
        when Symbol
          self.instance_eval <<-EOT
            def_var :#{blob} do |obj|
              (obj || job_item).property(:#{blob}) || unescape_dir_string("#{full_name}", obj)
            end
          EOT
        end
      end
    end
  end
  module BaseConfig
    extend Pipeline::Config

    def_var :next_step do @script.steps[ @script.steps.index(step) + 1 ] end
    def_var :prev_step do @script.steps[ @script.steps.index(step) - @script.steps.size - 1 ] end
    def_var :pipe do @script.class.name.sub(/::.*/,"").snake_case.to_sym end
    def_var :script do @script.class.name.snake_case.to_sym end
    def_var :pipe_script do "#{pipe}_#{script}" end


    def_var :logging_level do :WARN end

    def_var :error_file do "ERROR.#{cohort_name}" end
    def_var :verbose? do verbose == "yes" || verbose == "true" end
    def_var :config_dir do "#{lib_dir}/config" end
    def_var :tools_config do "#{config_dir}/tools.yml" end
    def_var :genome_config do "#{config_dir}/#{genome}.yml" end

    def_var :threads do resources[:threads] end
    def_var :qual_type do "phred64" end

    def_var :genome do :hg19 end
    def_var :reference_name do genome end
    def_var :reference_date do send "#{genome}_date".to_sym end
    def_var :reference_fa do send "#{genome}_fa".to_sym end
    def_var :reference_dict do send "#{genome}_dict".to_sym end
    def_var :reference_snp_vcf do send "#{genome}_snp_vcf".to_sym end
    def_var :reference_indel_vcf do send "#{genome}_indel_vcf".to_sym end
    def_var :reference_gtf do send "#{genome}_reference_gtf".to_sym end
    def_var :reference_2bit do send "#{genome}_2bit".to_sym end
    def_var :reference_rsem do send "#{genome}_rsem".to_sym end
    def_var :reference_interval_bed do send "#{genome}_interval_bed".to_sym end
    def_var :reference_gc do send "#{genome}_gc".to_sym end

    dir_tree({
      ":log_dir" => {
        ":pipe.:cohort_name.:step.:job_index.log" => :step_log,
        ":pipe.:cohort_name.:script.log" => :main_log
      },
      ":scratch_dir" => {
        "@sample_name" => {
          "." => :sample_scratch
        },
        "@cohort_name" => {
          "." => :cohort_scratch,
          "error.pid" => :error_pid
        }
      },
      ":metrics_dir" => {
        "timer_metrics" => :timer_metrics
      },
    })

    def_var :modules do [ :default ] end

    empty_var :work_dir, :verbose, :job_number, :keep_temp_files, :single_step, :filter_config, :walltime, :job_array

    def job_index sample=nil
      @opts[:job_number] ? @opts[:job_number] - 1 : 0 
    end
    def job_item sample=nil
      @opts[:job_array] ? @opts[:job_array][job_index] : @config
    end
    def trials
      job_array ? job_array.size : nil
    end
    def bunch unit
      send("#{unit}s".to_sym)
    end
    def find_job_item unit, item
      if item
        bunch(unit).find{ |unit| unit.property("#{unit}_name".to_sym) == item }
      else
        job_item.owner "#{unit}_name".to_sym
      end
    end


    def cohort
      return @config
    end

    job_items :sample, :replicate
    def_var :input_bam do |s| (s || job_item).property :input_bam end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do samples.map{ |s| sample_bam(s) } end
    def_var :sample_name do |s| (s || job_item).property :sample_name end
    def_var :sample_names do job_array.map{ |item| item.property :sample_name }.uniq end

    def_var :normal_bam do sample_bam(normal) end
    def_var :normal_name do |s| (s || job_item).property(:normal_name) || samples.first.sample_name end
    def_var :normal do samples.find{|s| s.sample_name == normal_name} end

    def normal_samples
      n = samples.map(&:normal_name).compact
      n.empty? ? [samples.first] : samples.select{|s| n.include? s.sample_name}
    end

    def tumor_samples
      n = samples.select { |s| s[:normal_name] }
      n.empty? ? samples[1..-1] : n
    end
    def_var :tumor_bam do sample_bam(sample) end

    def chromosomes
      # read chromosomes from the fasta dict
      @chroms ||= File.foreach(reference_dict).map do |s|
        s.match(/@SQ.*SN:(\w+)\s/) do |m| 
          next if m[1] == "chrM"
          { :chrom_name => m[1] }
        end
      end.compact
    end

    def initialize(script,opts=nil)
      @script = script
      # stuff that should default to nil goes here, rather than using method_missing
      @config = {}

      # also include the script config
      #self.class.send(:include, script.class.daughter_class(:config))

      load_opts
      load_procs
      load_env_vars
      load_tools_config

      set_dir

      load_config
      load_genome_config

      init_hook
    end

    def init_hook
    end

    def load_opts
      @opts = {}
      self.class.ancestors.each do |a|
        @opts.update a.class_opts if a.respond_to? :class_opts
      end
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

    def load_genome_config
      @opts.update( YAML.load_file( genome_config ) )
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

    def hostname
      @hostname ||= %x{ hostname }
    end

    def load_config
      @config=Pipeline::SampleObject.new YAML.load_file(config_file), self
    end

    def unescape_dir_string str, obj
      obj ||= job_item
      str.gsub(/:(?:(\w+)|\{(\w+)\})/) do |s|
        # send it with the appropriate argument
        send ($1 || $2).to_sym, obj
      end.gsub(/@(?:(\w+)|\{(\w+)\})/) do |s|
        m = ($1 || $2).to_sym
        obj.property m
      end
    end

    def find_prop property, item
      item ||= job_item
      item.property(property) || send(property, item)
    end

    def is_proc_collection? meth
      meth =~ /^\w+__\w+s$/
    end

    def proc_collect m, args
      terms = m.to_s.split(/__/)
      meth = terms.pop.sub!(/s$/,"").to_sym
      obj = [ job_item ]
      terms.each do |key|
        key = key.to_sym
        obj = obj.collect{|o| o.property key}.flatten
      end
      obj.map do |o|
        send meth, o
      end.flatten
    end

    def method_missing(meth,*args,&block)
      # always look in the config file first, so we can overrule things there
      meth = meth.to_sym if meth.is_a? String
      if @config[meth]
        return @config[meth]
      # check for it on the job_item
      elsif job_item && job_item != @config && job_item.property(meth)
        return job_item.property meth
      # if not there, it could be in the procs table, to be found interactively
      elsif @procs[meth]
        return instance_exec( *args, &@procs[meth])
      # if not there, it might be automatically collecting procs
      elsif is_proc_collection?(meth)
        return proc_collect(meth, args)
      # if not there, it is likely in options. These can be nil for empty variables
      elsif @opts.has_key? meth
        return @opts[meth]
      else
        super
      end
    end
  end
end
