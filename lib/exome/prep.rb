#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'

module Exome
  class Prep
    include Pipeline::Step
    runs_task :create_intervals_bed

    class CreateIntervalsBed
      include Pipeline::Task
      requires_file :interval_list
      dumps_file :interval_bed
      
      def run
        if !File.exists? config.interval_bed
          File.open(config.interval_bed,"w") do |f|
            File.foreach(config.interval_list) do |l|
              next if l =~ /^@/
              f.print l
            end
          end
        end
      end
    end
  end
end
