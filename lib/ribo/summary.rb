#!/usr/bin/env ruby

require 'flagstat'
require 'hash_table'
require 'gtf'

module Ribo
  class Summary
    include Pipeline::Step
    runs_tasks :summarize_normal_cov, :summarize_qc

    class SummarizeQc
      include Pipeline::Task
      outs_file :qc_summary
      
      def run
        qc = {}
        config.fractions.each do |f|
          # read in each file name and build up a hash of interesting information
          sample_qc = {}
          flags = Flagstat.new config.qc_flag(f)
          flags.each do |flag,value|
            sample_qc[flag] = value.first
          end

          sample_qc[:splice_counts] = File.read(config.qc_splice_counts f).to_i
          qc[s.sample_name] = sample_qc

          align = HashTable.new comment: /^(#|$)/
          align.parse config.qc_align_metrics(f)
          align.columns.each do |flag|
            sample_qc[flag] = align.first[flag]
          end
          
          rnaseq = HashTable.new comment: /^(#|$)/, downcase: true
          rnaseq.parse config.qc_rnaseq(s)
          rnaseq.columns.each do |flag|
            sample_qc[flag] = rnaseq.first[flag]
          end
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
  end
end
