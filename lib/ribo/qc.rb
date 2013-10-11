#!/usr/bin/env ruby
module Ribo
  class Qc
    include Pipeline::Step
    runs_tasks :calc_flags, :calc_rna_metrics, :collect_align_metrics, :collect_splice_counts
    runs_on :samples

    class CalcFlags
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :qc_flag

      def run
	log_info "Calculate flag statistics"
        sam_flags config.sample_bam, config.qc_flag or error_exit "Flagstats failed"
      end
    end

    class CalcRnaMetrics
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :qc_rnaseq, :qc_pdf

      def run
	log_info "Calculate rnaseq metrics"
        picard :collect_rna_seq_metrics, :ASSUME_SORTED => :true, :REF_FLAT => config.hg19_refflat, :RIBOSOMAL_INTERVALS => config.hg19_rrna_intervals, :STRAND_SPECIFICITY => :NONE,
          :CHART_OUTPUT => config.qc_pdf, :INPUT => config.sample_bam, :OUTPUT => config.qc_rnaseq or error_exit "RNAseq metrics failed"
      end
    end

    class CollectAlignMetrics
      include Pipeline::Task
      requires_file :sample_bam
      outs_files :qc_align_metrics

      def run
        log_info "Calculating alignment metrics"
        picard :collect_alignment_summary_metrics, :INPUT => config.sample_bam, :OUTPUT => config.qc_align_metrics or error_exit "Alignment metrics failed"
      end
    end

    class CollectSpliceCounts
      include Pipeline::Task
      requires_file :sample_bam
      outs_files :qc_splice_counts

      def run
        run_cmd %Q!samtools view #{config.sample_bam} | awk -F '	' '$6 ~ /N/ { print "OK"}' | wc -l > #{config.qc_splice_counts}! or error_exit "Couldn't count spliced reads"
      end
    end
  end
end
