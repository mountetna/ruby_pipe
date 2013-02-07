require 'readline'

module Pipeline
  module ConfigGenerator
    # Class methods and variables
    module ClassMethods
      attr_reader :usages

      def usage cmd, expln
        @usages ||= {
          :reset => [ nil, 0, "create a fresh config file" ],
          :show =>  [ nil, 0, "list the current config file" ],
          :print => [ nil, 0, "print the current config to a local file and exit" ],
          :samples => [ "<list of sample names>", 1, "set the sample names" ]
        }
        cmd,args = cmd.split(/ /,2)
        @usages.update Hash[ cmd.to_sym, [ args, required_args(args), expln ] ]
      end

      def required_args cmd
        cmd.gsub(/\[.*?\]/,"").scan(/<.*?>/).size
      end
    end

    def self.included(base)
      base.extend ClassMethods
    end

    # Instance methods and variables
    attr_reader :config

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

    def default_file
      "#{@config[:job_name]}.#{parent_class}.yml"
    end

    def commands
      [ :samples, :show, :print, :reset ] + public_methods(nil)
    end

    def initialize j
      if j =~ /\.(yaml|yml|conf)$/
        @config_file = j
      else
        @job_name = j
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
      while l = Readline.readline(":".red,false)
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

    def usage
      puts "Commands:"
      self.class.usages.each do |c,u|
        cmd = [ c, u.first ].compact.join " "
        puts " %-50s" % cmd.bold + "# #{u.last}".cyan
      end
      nil
    end

    def reset(args=nil)
      @config = { :job_name => @job_name}.update YAML.load( config_text )
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

    def print args=nil
      File.open(args.first || default_file, "w") do |f|
        f.puts config_header
        f.puts
        f.puts config.to_yaml
        f.puts
        f.puts config_comments
      end
    end

    def check_usage cmd, args
      required = self.class.usages[cmd][1]
      if args.size < required
        usage
        return true
      end
    end

    def start_generator
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
