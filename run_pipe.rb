#!/usr/bin/env ruby

ENV['LIB_DIR'] = '/taylorlab/scripts/ruby_pipe'
  #File.dirname(File.realdirpath(__FILE__))

$: << ENV['LIB_DIR'] + "/lib"

pipe,script = File.basename(__FILE__).split(/\_/,2)

pipe.sub! /^local./,""

require pipe + '/' + script

s = Kernel.const_get( pipe.camel_case).const_get( script.camel_case).new

begin
  s.init ARGV
rescue Errno::EPIPE
  exit
end
