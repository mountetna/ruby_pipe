#!/usr/bin/env ruby
module Genome
  class CollectQc
  #class Qc
    include Pipeline::Step
    runs_tasks :calc_flags, :collect_insert_sizes, :collect_align_metrics #, :coverage_metrics
    runs_on :samples

    class CalcFlags
      include Pipeline::Task
      requires_file :qc_bam
      outs_file :qc_flag

      def run
	log_info "Calculate flag statistics"
        sam_flags config.qc_bam, config.qc_flag or error_exit "Collecting flag statistics failed"
      end
    end

    class CollectInsertSizes
      include Pipeline::Task
      requires_file :qc_bam
      outs_files :qc_inserts, :qc_histogram

      def run
        log_info "Calculating insert metrics"
        picard :collect_insert_size_metrics, :INPUT => config.qc_bam, :OUTPUT => config.qc_inserts, :HISTOGRAM_FILE => config.qc_histogram or error_exit "Collecting insert sizes failed"
      end
    end

    class CollectAlignMetrics
      include Pipeline::Task
      requires_file :qc_bam
      outs_files :qc_align_metrics

      def run
        log_info "Calculating alignment metrics"
        picard :collect_alignment_summary_metrics, :INPUT => config.qc_bam, :OUTPUT => config.qc_align_metrics or error_exit "Collecting alignment metrics failed"
      end
    end

    class CoverageMetrics
      include Pipeline::Task
      requires_file :qc_bam
      outs_file :qc_coverage_metrics

      def run
        #create_interval_bed
        gatk :depth_of_coverage, :out => config.qc_coverage_base, :input_file =>  config.qc_bam,
          :omitDepthOutputAtEachBase => true,
          :omitIntervalStatistics => true,
          :omitLocusTable => true
          #:intervals => config.interval_bed or error_exit "Coverage computation failed"
      end
    end

  end

  #class QcSummary
  class CollectQcSummary
    include Pipeline::Step
    runs_tasks :summarize_qc

    class SummarizeQc
      include Pipeline::Task
      outs_file :qc_summary
      
      def run
        qc = {}
        config.samples.each do |s|
          # read in each file name and build up a hash of interesting information
          sample_qc = {}
          flags = Hash[[ :total, :duplicates, :mapped, :paired_in_sequence, :read1, :read2, :properly_paired, :both_mapped, :singletons, :mate_mapped_chr, :mate_mapped_chr_highq ].zip File.foreach(config.qc_flag s).map(&:split).map(&:first)]
          sample_qc.update flags

          #hybrid = File.foreach(config.qc_hybrid s).reject{|i| i=~ /^(#|$)/}.map(&:split)
          #hybrid = Hash[hybrid.first.map(&:downcase).map(&:to_sym).zip hybrid.last]
          #sample_qc.update hybrid

          #config.qc_align_metrics - these are maybe redundant?
          
          inserts = File.foreach(config.qc_inserts s).reject{|i| i=~ /^(#|$)/}.map(&:split)[0..1]
          inserts = Hash[inserts.first.map(&:downcase).map(&:to_sym).zip inserts.last]
          sample_qc[:median_insert_size] = inserts[:median_insert_size]
          sample_qc[:mean_insert_size] = inserts[:mean_insert_size]

          #coverage = File.foreach(config.qc_coverage_metrics s).reject{|i| i=~ /^(#|$)/}.map(&:split)[0..1]
          #coverage = Hash[coverage.first.map(&:downcase).map(&:to_sym).zip coverage.last]
          #sample_qc[:mean_coverage] = coverage[:mean]

          qc[s.sample_name] = sample_qc
        end
        if qc.size > 1
	    log_info "Inside if qc.save 1"
            File.open(config.qc_summary, "w") do |f|
            metrics = qc.first.last.keys
            samples = qc.keys
            f.puts "\t" + samples.join("\t")
            metrics.each do |metric|
              f.puts "#{metric}\t#{samples.map{|s| qc[s][metric]}.join("\t")}"
            end
          end
        end
      end
    end
  end
end
