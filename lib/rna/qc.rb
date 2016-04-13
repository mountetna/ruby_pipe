require 'picard_metrics'
require 'flagstat'

module Rna
  class Qc
    include Pipeline::Step
    runs_tasks :calc_flags, :calc_rna_metrics, :collect_align_metrics
    runs_on :samples, :replicates

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
  class QcSummary
    include Pipeline::Step
    runs_on :cohort
    runs_tasks :summarize_qc

    class SummarizeQc
      include Pipeline::Task
      requires_files :samples__replicates__qc_rnaseqs
      requires_files :samples__replicates__qc_flags
      outs_file :qc_summary

      def run
        table = nil
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            mets = PicardMetrics.new
            mets.parse config.qc_rnaseq(rep)
            flags = Flagstat.new config.qc_flag(rep)

            if !table
              table = HashTable.new columns: [:replicate] + 
                flags.clean_hash.keys + 
                mets.sections[:rna_seq_metrics].metrics.keys
            end
            table << mets.sections[:rna_seq_metrics].metrics.merge(
              replicate: config.sample_replicate_name(rep)
            ).merge( flags.clean_hash )
          end
        end
        table.print config.qc_summary
      end
    end
  end
end
