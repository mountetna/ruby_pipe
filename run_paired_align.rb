#!/usr/bin/env ruby

ENV['LIB_DIR'] = File.dirname(__FILE__)
$: << ENV['LIB_DIR'] + "/lib"

# make some dummy variables

require 'exome/paired_align'

# Create the pipe
pipe = Exome::PairedAlign.new

# do what needs to be done
pipe.init ARGV
