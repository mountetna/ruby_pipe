#!/usr/bin/env ruby
module Exome
  class Merge
    include Pipeline::Step
    runs_tasks :merge_sam_files, :index_bam
    job_list do config.samples end
    resources :threads => 12

    class MergeSamFiles
      include Pipeline::Task
      requires_file :aligned_bams
      outs_file :raw_sample_bam

      def run
        log_info "Merging sample bam files for #{config.sample_name}"
        picard :merge_sam_files, :CREATE_INDEX => :true, :O => config.raw_sample_bam, :I => config.aligned_bams or error_exit "BAM file merge failed."
      end
    end
    class IndexBam
      include Pipeline::Task
      requires_file :raw_sample_bam
      outs_file :raw_sample_bai

      def run
        log_info "Indexing sample bam file."
        sam_index config.raw_sample_bam or error_exit "Bam indexing failed."
      end
    end
  end
end
