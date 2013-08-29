#!/usr/bin/env ruby
#
require 'yaml'
require 'readline'

exit if !ARGV[0] || !File.exists?(ARGV[0])

ENV['LIB_DIR'] = File.dirname File.realdirpath(__FILE__)
ENV['CONFIG'] = ARGV[0]
__config__ = YAML.load_file ENV['CONFIG']
$: << ENV['LIB_DIR'] + "/lib"

# make some dummy variables

require __config__[:pipe].downcase + '/' + __config__[:script]
# Create the pipe
pipe = Kernel.const_get( __config__[:pipe].camel_case ).const_get( __config__[:script].camel_case).new

# do what needs to be done
ARGV.clear

class SampleGenerator
  include Pipeline::Usage
  PIPELINE_SAMPLE_FILE = "/taylorlab/data/timur_data/pipeline_samples.yml"

  def commands
    usages.keys
  end

  def pipeline_samples
    @pipeline_samples ||= YAML.load_file(PIPELINE_SAMPLE_FILE)|| []
  end

  def show args=nil
    puts pipeline_samples.to_yaml
  end
  usage "show", "list the current config file"

  def setup_generator
    Readline.completion_proc = Proc.new do |str|
        opts = Dir[str+'*'].grep( /^#{Regexp.escape(str)}/ )
        if opts.length == 1 && File.directory?(opts.first)
          Readline.completion_append_character = "/"
        else
          Readline.completion_append_character = " "
        end
        opts
    end
  end

  class Sample
    def get_line prop
      Readline.readline("#{prop}: ".red.bold,true)
    end
    def initialize barcode
      @sample = {}
      @sample[:barcode] = barcode
    end
    def get_prop prop
      @sample[prop] = get_line prop
    end
    def encode_with coder
      coder.tag = nil
      @sample.each do |k,v|
        coder[k.to_s] = v
      end
    end

    def valid?
      @sample[:description] && @sample[:lab] && @sample[:sample_type] && @sample[:cancer_type]
    end
  end

  def start_generator
    while buf = Readline.readline("> ",true)
      exit if !buf
      cmd, buf = buf.split
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

  def list_samples
    puts "Found samples:".blue.bold, @samples
  end

  def update_sample sample_name
    return nil if !sample_name
    puts "Updating #{sample_name}".blue.bold
    sample = Sample.new sample_name
    sample.get_prop :description or return nil
    sample.get_prop :lab or return nil
    sample.get_prop :sample_type or return nil
    sample.get_prop :cancer_type or return nil
    if sample.valid?
      pipeline_samples.push sample
    end
  end

  def update args=nil
    if !args
      @samples.shift if update_sample @samples.first
    else
      update_sample @samples.find{|i| i == args.first}
    end
  end
  usage "update [<sample>]", "Add information for sample or first name"

  def prune_existing_samples
    @samples.reject! do |s|
      pipeline_samples.find{|e| e['barcode'] == s }
    end
  end

  def print args=nil
    File.open PIPELINE_SAMPLE_FILE, "w" do |f|
      f.puts pipeline_samples.to_yaml
    end
  end
  usage "print", "print the current config file to a local file"

  def initialize ss
    @samples = ss

    prune_existing_samples

    list_samples

    setup_generator
  end
end

gen = SampleGenerator.new pipe.config.samples.map(&:timur_name).compact
gen.start_generator
