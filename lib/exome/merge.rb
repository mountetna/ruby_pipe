#!/usr/bin/env ruby
module Exome
  class Merge
    include Pipeline::Step
    runs_tasks :merge_sam_files
    runs_on :samples
    resources :threads => 12

    class MergeSamFiles
      include Pipeline::Task
      requires_file :aligned_bams
      outs_file :raw_sample_bam

      def run
        log_info "Merging sample bam files for #{config.sample_name}"
        picard :merge_sam_files, :O => config.raw_sample_bam, :I => config.aligned_bams or error_exit "BAM file merge failed."
      end
    end
  end
end
