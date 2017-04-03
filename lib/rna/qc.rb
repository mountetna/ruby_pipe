require 'picard_metrics'
require 'flagstat'

module Rna
  class Qc
    include Pipeline::Step
    runs_tasks :calc_flags, :calc_rna_metrics, :collect_align_metrics, :bwa_calc_flags, :bwa_calc_rna_metrics #, :compute_exon_duplicate_rate
    runs_on :samples, :replicates
    resources memory: "20gb"

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
        picard :collect_rna_seq_metrics, :ASSUME_SORTED => :true, :REF_FLAT => config.reference_refflat, :RIBOSOMAL_INTERVALS => config.reference_rrna_intervals, :STRAND_SPECIFICITY => :NONE,
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

    class BwaCalcFlags
      include Pipeline::Task
      requires_file :bwa_aligned_bam
      outs_file :bwa_qc_flag

      def run
	log_info "Calculate flag statistics"
        sam_flags config.bwa_aligned_bam, config.bwa_qc_flag or error_exit "Flagstats failed"
      end
    end

    class BwaCalcRnaMetrics
      include Pipeline::Task
      requires_file :bwa_aligned_bam
      outs_file :bwa_qc_rnaseq, :bwa_qc_pdf

      def run
	log_info "Calculate rnaseq metrics"
        picard :collect_rna_seq_metrics, :ASSUME_SORTED => :true, :REF_FLAT => config.reference_refflat, :RIBOSOMAL_INTERVALS => config.reference_rrna_bwa_intervals, :STRAND_SPECIFICITY => :NONE,
          :CHART_OUTPUT => config.bwa_qc_pdf, :INPUT => config.bwa_aligned_bam, :OUTPUT => config.bwa_qc_rnaseq or error_exit "RNAseq metrics failed"
      end
    end

    class ComputeExonDuplicateRate
      include Pipeline::Task
      requires_file :output_bam
      outs_files :exon_unique_cov, :exon_total_cov

      def run
        log_info "Mapping coverage to reference genes"

        run_cmd "samtools view -h #{config.output_bam} > #{config.total_sam}"

        htseq_count input: config.total_sam, 
          gtf: config.flat_reference_gtf, 
          type: "exon",
          idattr: "exon_id",
          order: "pos",
          out: config.exon_total_cov or error_exit "Computing coverage for total failed."
        File.unlink config.total_sam

        run_cmd "samtools view -h -F 1024 #{config.output_bam} > #{config.unique_sam}"

        htseq_count input: config.unique_sam, 
          gtf: config.flat_reference_gtf, 
          type: "exon",
          idattr: "exon_id",
          order: "pos",
          out: config.exon_unique_cov or error_exit "Computing coverage for unique failed."
        File.unlink config.unique_sam
      end
    end
  end
 
  
  class BwaQcSummary
    include Pipeline::Step
    runs_on :cohort
    runs_tasks :bwa_summarize_qc
    resources memory: "10gb"
    
    class SummarizeQc
      include Pipeline::Task
      requires_files :samples__replicates__bwa_qc_rnaseqs
      requires_files :samples__replicates__bwa_qc_flags
      outs_file :bw_qc_summary

      def run
        table = nil
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            mets = PicardMetrics.new
            mets.parse config.bwa_qc_rnaseq(rep)
            flags = Flagstat.new config.bwa_qc_flag(rep)

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
        table.print config.bwa_qc_summary
      end
    end
  end

  class QcSummary
    include Pipeline::Step
    runs_on :cohort
    runs_tasks :summarize_qc
    resources memory: "10gb"
    
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
