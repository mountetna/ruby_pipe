#!/usr/bin/env ruby

ENV['ACTION'] = "exec"
ENV['STEP'] = "align"
ENV['LIB_DIR'] = File.dirname(__FILE__)

$: << "#{ENV['LIB_DIR']}/lib"

# make some dummy variables

require 'exome/paired_align'

# Create the current step
pipe = Exome::PairedAlign.new

# do what needs to be done
pipe.run_action
