#!/usr/bin/env ruby
require 'tempfile'

module Rna
  class AssembleTranscripts
    include Pipeline::Step
    runs_tasks :compare, :merge
    resources :threads => 12
    job_list do [ config.cohort ] end

    class Compare
      include Pipeline::Task
      requires_file :transcripts_gtfs
      dumps_file :tracking_file

      def run
        File.open( config.assembly_list, "w" ) do |f|
          f.puts config.transcripts_gtfs.join("\n")
        end

        log_info "Comparing transcripts to reference"
        cuffcompare :i => config.assembly_list, :r => config.hg19_ucsc_gtf, :o => config.cuffcompare_scratch

        File.unlink config.assembly_list
      end
    end
    class Merge
      include Pipeline::Task
      requires_files :transcripts_gtfs
      outs_file :assembly_gtf

      def run
        File.open( config.assembly_list, "w" ) do |f|
          f.puts config.transcripts_gtfs.join("\n")
        end

        log_info "Merging transcripts"
        cuffmerge :out => config.scratch, :list => config.assembly_list or error_exit "Cuffmerge failed."

        File.unlink config.assembly_list
      end
    end
  end
end
