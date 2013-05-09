require 'readline'

module Pipeline
  module ConfigGenerator
    include Pipeline::Usage
    attr_reader :config

    usage "reset", "create a fresh config file"
    usage "show", "list the current config file"
    usage "print", "print the currante config file to a local file"
    usage "samples <list of sample names>", "set the sample names"
    usage "quit", "quit the config generator"

    def lib_dir
      ENV['LIB_DIR']
    end

    def parent_class
      self.class.parent_class.class_symbol
    end

    def config_file
      @config_file || default_config
    end

    def config_text
      @config_text ||= File.read config_file
    end

    def default_config
      File.join lib_dir, "config", "#{parent_class}.yml"
    end

    def output_file
      @config_file || default_output
    end

    def default_output
      "#{@config[:cohort_name]}.#{parent_class}.yml"
    end

    def commands
      usages.keys
    end

    def initialize j
      if j =~ /\.(yaml|yml|conf)$/
        @config_file = j
      else
        @cohort_name = j
      end

      reset

      start_generator
    end

    def pair_fastqs fastq_list
      # for each fastq, find the pair in this list
      fastq_list.group_by do |s|
        s.sub(/_R[12][_.]/,"") 
      end.values.map do |s|
        Hash[[ :fq1, :fq2 ].zip s] 
      end
    end

    def comment txt
      puts txt.cyan
    end

    def getlines txt
      comment "List #{txt} separated by space or newline."
      comment "End with a . on a line by itself or ctrl-d to finish."
      while l = Readline.readline(":".red,true)
        break if l == "."
        if l =~ /^\!/
          cmd = l.sub(/^\!/,"")
          Readline::HISTORY << l
          puts out = `#{cmd}`
          yield out
        else
          yield l
        end
      end
      puts
    end

    def reset(args=nil)
      @config = { :cohort_name => @cohort_name}.update(YAML.load( config_text ) || {})
    end

    def samples(args=nil)
      config[:samples] = []
      args.each do |s|
        config[:samples].push({ :sample_name => s })
      end
    end

    def find_sample sample_name
      config[:samples] ||= []
      config[:samples].find { |s| s[:sample_name] == sample_name } 
    end

    def get_sample sample_name
      sample = find_sample sample_name
      if !sample
        sample = { :sample_name => sample_name }
        config[:samples].push sample
      end
      return sample
    end

    def config_header
      config_text.split(/\n/).grep(/^\s*#/).first.chomp
    end

    def config_comments
      config_text.split(/\n/).grep(/^\s*#/)[1..-1].map(&:chomp)
    end

    def show args=nil
      puts config_header.green
      puts
      puts config.to_yaml
      puts
      puts config_comments.map(&:green)
    end

    def quit args=nil
      exit
    end

    def print args=nil
      File.open(args.first || default_output, "w") do |f|
        f.puts config_header
        f.puts
        f.puts config.to_yaml
        f.puts
        f.puts config_comments
      end
    end

    def start_generator
      Readline.completion_proc = Proc.new do |str|
          opts = Dir[str+'*'].grep( /^#{Regexp.escape(str)}/ )
          if opts.length == 1 && File.directory?(opts.first)
            Readline.completion_append_character = "/"
          else
            Readline.completion_append_character = " "
          end
          opts
      end
      while buf = Readline.readline("> ",true)
        exit if !buf
        buf = buf.split
        cmd = buf.shift
        next if !cmd
        cmd = cmd.to_sym
        if !commands.include? cmd
          usage
          next
        end
        next if check_usage cmd, buf
        send cmd, buf
      end
      puts
    end
  end
end
