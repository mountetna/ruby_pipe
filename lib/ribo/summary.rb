#!/usr/bin/env ruby

require 'flagstat'
require 'hash_table'
require 'gtf'

module Ribo
  class Summary
    include Pipeline::Step
    runs_tasks :summarize_qc, :summarize_normal_cov, :summarize_null_cov

    class SummarizeQc
      include Pipeline::Task
      outs_file :qc_summary
      
      def run
        qc = {}
        config.samples.each do |s|
          # read in each file name and build up a hash of interesting information
          sample_qc = {}
          flags = Flagstat.new config.qc_flag(s)
          flags.each do |flag,value|
            sample_qc[flag] = value.first
          end

          sample_qc[:splice_counts] = File.read(config.qc_splice_counts s).to_i
          qc[s.sample_name] = sample_qc

          align = HashTable.new config.qc_align_metrics(s), :comment => /^(#|$)/, :downcase => true
          align.header.each do |flag|
            sample_qc[flag] = align.first[flag]
          end
          
          rnaseq = HashTable.new config.qc_rnaseq(s), :comment => /^(#|$)/, :downcase => true
          rnaseq.header.each do |flag|
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
      
      def run
        summary = {}
        gtf = GTF.new config.hg19_unified_gtf, :idx => [ :gene_id ]
        config.samples.each do |s|
          # read in each file name and build up a hash of interesting information
          cov = HashTable.new config.normal_cov(s), :header => [ :gid, :count ]
          cov.each do |g|
            summary[g.gid] ||= {}
            summary[g.gid][s.sample_name] = g.count
          end
        end
        File.open(config.normal_summary,"w") do |f|
          f.puts [ :gene_id, :symbol, :size, config.samples.map(&:sample_name) ].flatten.join("\t")
          summary.each do |gid,samples|
            next if gid !~ /ENS/
            symbol = gtf.idx(:gene_id,gid).first.gene_name
            gene_size = gtf.idx(:gene_id,gid).inject(0) {|sum,f| sum += f.size }
            f.puts [ gid, symbol, gene_size, config.samples.map{|s| samples[s.sample_name]} ].flatten.join("\t")
          end
        end
      end
    end
    class SummarizeNullCov
      include Pipeline::Task
      outs_file :null_summary

      def run
        summary = {}
        config.samples.each do |s|
          # read in each file name and build up a hash of interesting information
          cov = HashTable.new config.null_cov(s), :header => [ :chr, :type, :feature, :start, :stop, :score, :strand, :frame, :att, :count ]
          cov.each do |g|
            key = [ g.chr, g.start, g.stop ]
            summary[key] ||= {}
            summary[key][s.sample_name] = g.count
          end
        end
        File.open(config.null_summary,"w") do |f|
          f.puts [ :gene_id, :symbol, :size, config.samples.map(&:sample_name) ].flatten.join("\t")
          summary.each do |key,samples|
            f.puts [ "#{key[0]}:#{key[1]}-#{key[2]}", :random, key[2].to_i-key[1].to_i, config.samples.map{|s| samples[s.sample_name]} ].flatten.join("\t")
          end
        end
      end
    end
  end
end
