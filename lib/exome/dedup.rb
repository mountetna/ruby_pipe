#!/usr/bin/env ruby
module Exome
  class Dedup
    include Pipeline::Step
    runs_tasks :mark_duplicates
    runs_on :samples
    resources :threads => 12

    class MarkDuplicates
      include Pipeline::Task
      requires_file :aligned_bams
      outs_file :raw_sample_bam, :duplication_metrics

      def run
	log_info "Mark duplicates"
	picard :mark_duplicates, :INPUT => config.aligned_bams,
          :OUTPUT => config.raw_sample_bam,
          :METRICS_FILE => config.duplication_metrics, 
          :CREATE_INDEX => :true,
          :REMOVE_DUPLICATES => :false or error_exit "Mark duplicates failed"
      end
    end
  end
end
