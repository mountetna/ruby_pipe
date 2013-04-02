#!/usr/bin/env ruby
module Rna
  class Qc
    include Pipeline::Step
    runs_tasks :calc_flags, :calc_rna_metrics, :collect_align_metrics
    job_list do config.replicates end

    class CalcFlags
      include Pipeline::Task
      requires_file :qc_bam
      outs_file :qc_flag

      def run
	log_info "Calculate flag statistics"
        sam_flags config.qc_bam, config.qc_flag or error_exit "Flagstats failed"
      end
    end

    class CalcRnaMetrics
      include Pipeline::Task
      requires_file :qc_bam
      outs_file :qc_rnaseq, :qc_pdf

      def run
	log_info "Calculate rnaseq metrics"
        picard :collect_rna_seq_metrics, :ASSUME_SORTED => :true, :REF_FLAT => config.hg19_refflat, :RIBOSOMAL_INTERVALS => config.hg19_rrna_intervals, :STRAND_SPECIFICITY => :NONE,
          :CHART_OUTPUT => config.qc_pdf, :INPUT => config.qc_bam, :OUTPUT => config.qc_rnaseq or error_exit "RNAseq metrics failed"
      end
    end

    class CollectAlignMetrics
      include Pipeline::Task
      requires_file :qc_bam
      outs_files :qc_align_metrics

      def run
        log_info "Calculating alignment metrics"
        picard :collect_alignment_summary_metrics, :INPUT => config.qc_bam, :OUTPUT => config.qc_align_metrics or error_exit "Alignment metrics failed"
      end
    end
  end
end
