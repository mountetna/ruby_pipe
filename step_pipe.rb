#!/usr/bin/env ruby

$: << ENV['LIB_DIR'] + "/lib"

# make some dummy variables

require ENV['PIPE'].downcase + '/' + ENV['SCRIPT']

# Create the pipe
pipe = Kernel.const_get( ENV['PIPE'].camel_case ).const_get( ENV['SCRIPT'].camel_case).new

# do what needs to be done
pipe.run_action
