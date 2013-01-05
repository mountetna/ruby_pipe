#!/usr/bin/env ruby
module Exome
  class HybridQc
    include Pipeline::Step
    runs_tasks :calc_flags, :calc_metrics, :collect_insert_sizes, :collect_align_metrics

    class CalcFlags
      include Pipeline::Task
      requires_file :qc_bam
      outs_file :qc_flag

      def run
	log_info "Calculate flag statistics"
        sam_flags config.qc_bam, config.qc_flag
      end
    end

    class CalcMetrics
      include Pipeline::Task
      requires_file :qc_bam
      outs_file :qc_hybrid

      def run
	log_info "Calculate hybrid selection metrics"
        picard :calculate_hs_metrics, :BAIT_INTERVALS => config.interval_list, :TARGET_INTERVALS => config.interval_list, :INPUT => config.qc_bam, :OUTPUT => config.qc_hybrid
      end
    end

    class CollectInsertSizes
      include Pipeline::Task
      requires_file :qc_bam
      outs_files :qc_inserts, :qc_histogram

      def run
        log_info "Calculating insert metrics"
        picard :collect_insert_size_metrics, :INPUT => config.qc_bam, :OUTPUT => config.qc_inserts, :HISTOGRAM_FILE => config.qc_histogram
      end
    end

    class CollectAlignMetrics
      include Pipeline::Task
      requires_file :qc_bam
      outs_files :qc_align_metrics

      def run
        log_info "Calculating alignment metrics"
        picard :collect_alignment_summary_metrics, :INPUT => config.qc_bam, :OUTPUT => config.qc_align_metrics
      end
    end
  end
end
