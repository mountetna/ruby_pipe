#!/usr/bin/env ruby
module Exome
  class Merge
    include Pipeline::Step
    runs_tasks :merge_sam_files, :index_bam
    class MergeSamFiles
      include Pipeline::Task
      requires_file :aligned_bams
      outs_file :merged_bam

      def run
        log_info "Merging sample bam files for #{config.sample_name}"
        picard :merge_sam_files, :CREATE_INDEX => :true, :O => config.merged_bam, :I => config.merged_bam or error_exit "BAM file merge failed."
      end
    end
    class IndexBam
      include Pipeline::Task
      requires_file :merged_bam
      outs_file :merged_bai

      def run
        log_info "Indexing sample bam file."
        sam_index config.merged_bam or error_exit "Bam indexing failed."
      end
    end
  end
end
