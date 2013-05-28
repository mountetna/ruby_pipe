#!/usr/bin/env ruby
#
require 'pry'
require 'yaml'

ENV['LIB_DIR'] = File.dirname File.realdirpath(__FILE__)
puts ENV['LIB_DIR']
ENV['CONFIG'] = ARGV[0]

__config__ = YAML.load_file ENV['CONFIG']

$: << ENV['LIB_DIR'] + "/lib"

# make some dummy variables

require __config__[:pipe].downcase + '/' + __config__[:script]
# Create the pipe
pipe = Kernel.const_get( __config__[:pipe].camel_case ).const_get( __config__[:script].camel_case).new


# do what needs to be done
ARGV.clear
binding.pry

