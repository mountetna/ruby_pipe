#!/usr/bin/env ruby

require 'picard_metrics'
require 'flagstat'
require 'hash_table'
require 'gtf'

module Ribo
  class Summary
    include Pipeline::Step
    runs_tasks :summarize_normal_cov, :summarize_transcript_model_cov, :summarize_qc

    class SummarizeQc
      include Pipeline::Task
      outs_file :qc_summary
      
      def run
        qc = {}
        config.fractions.each do |f|
          # read in each file name and build up a hash of interesting information
          sample_qc = {}
          qc[f.fraction_name] = sample_qc

          flags = Flagstat.new config.qc_flag(f)
          flags.each do |flag,value|
            sample_qc[flag] = value.first
          end

          mets = PicardMetrics.new
          mets.parse config.qc_rnaseq(f)

          sample_qc.update mets.sections[:rna_seq_metrics].metrics

          rrna_metrics = HashTable.new columns: [ :name, :count ], 
                          types: { count: :int }, 
                          index: [ :name ], 
                          parse_mode: :noheader
          rrna_metrics.parse config.qc_rrna_metrics f
          sample_qc[:rrna_reads] = rrna_metrics.index[:name]["rRNA_reads"].first.count
          sample_qc[:non_rrna_reads] = rrna_metrics.index[:name]["non_rRNA_reads"].first.count
          sample_qc[:pct_rrna] = (sample_qc[:rrna_reads] / (sample_qc[:rrna_reads] + sample_qc[:non_rrna_reads]).to_f).round(5)
        end
        if qc.size > 1
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
    class SummarizeNormalCov
      include Pipeline::Task
      outs_file :normal_summary

      class CoverageTable <  HashTable
        columns :gid, :count
        parse_mode :noheader
      end
      
      def run
        gtf = GTF.new index: [ :gene_id ]
        gtf.parse config.reference_unified_gtf
        summary = HashTable.new columns: [ :gene_id, :symbol, :size, config.fractions.map{|f|f.fraction_name.to_sym} ].flatten, index: [ :gene_id ]
        config.fractions.each do |s|
          # read in each file name and build up a hash of interesting information
          cov = CoverageTable.new.parse config.normal_cov(s)
          cov.each do |gene|
            if summary.index[:gene_id][gene.gid].count == 0
              summary << { 
                gene_id: gene.gid, 
                symbol: gtf.index[:gene_id][gene.gid].map(&:gene_name).first,
                size: gtf.index[:gene_id][gene.gid].inject(0) {|sum,f| sum += f.size },
                s.fraction_name.to_sym => gene.count
              }
            else
              entry = summary.index[:gene_id][gene.gid].first
              entry.update s.fraction_name.to_sym => gene.count
            end
          end
        end

        summary.print config.normal_summary
      end
    end
    class SummarizeTranscriptModelCov
      include Pipeline::Task
      requires_file :transcript_model_gtf
      outs_file :transcript_model_summaries

      class CoverageTable <  HashTable
        columns :gid, :count
        parse_mode :noheader
      end
      
      def run
        log_info "Loading GTF"
        gtf = GTF.new index: [ :gene_id ]
        gtf.parse config.transcript_model_gtf

        summary_columns = [ :gene_id, :symbol, :size, config.fractions.map{|f|f.fraction_name.to_sym} ].flatten
        
        config.transcript_model_regions.each do |region|
          summary = HashTable.new columns: summary_columns, index: [ :gene_id ]
          config.fractions.each do |s|
            # read in each file name and build up a hash of interesting information
            cov = CoverageTable.new.parse config.transcript_model_coverage(region, s)
            cov.each do |gene|
              if summary.index[:gene_id][gene.gid].count == 0
                summary << { 
                  gene_id: gene.gid, 
                  symbol: gtf.index[:gene_id][gene.gid].map(&:gene_name).first,
                  size: gtf.index[:gene_id][gene.gid].select{|f| f.feature == region.to_s}.inject(0){|sum,f| sum += f.size },
                  s.fraction_name.to_sym => gene.count
                }
              else
                entry = summary.index[:gene_id][gene.gid].first
                entry.update s.fraction_name.to_sym => gene.count
              end
            end
          end
          summary.print config.transcript_model_summary(region)
        end
      end
    end
  end
end
