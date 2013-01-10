#!/usr/bin/env ruby
require 'hash_table'
module Exome
  class CopyNumber
    include Pipeline::Step
    runs_tasks :compute_coverage, :compute_ratio
    class ComputeCoverage
      include Pipeline::Task
      requires_files :tumor_bam, :normal_bam, :interval_bed
      dumps_files :tumor_cov, :normal_cov

      def run
        coverage_bed config.normal_bam, config.interval_bed, config.normal_cov or error_exit "Computing normal coverage failed."
        coverage_bed config.tumor_bam, config.interval_bed, config.tumor_cov or error_exit "Computing tumor coverage failed."
      end
    end
    class ComputeRatio
      include Pipeline::Task
      requires_files :normal_cov, :tumor_cov, :insert_mutations
      outs_files :tumor_gene_cnr, :tumor_exon_cnr, :tumor_mutations

      def exon_coverage(cov)
        tot = cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        cov.each_with_object(Hash.new(0)) do |l,m| 
          m[l[:name]] += (l[:count].to_i + 10) / tot 
        end
      end

      def run
        normal_cov = HashTable.new(config.normal_cov, [ :chr, :start, :stop, :strand, :name, :count ])
        tumor_cov = HashTable.new(config.tumor_cov, [ :chr, :start, :stop, :strand, :name, :count ])
        tumor_logr = HashTable.new(config.tumor_cov, [ :chr, :start, :stop, :strand, :name, :count ])

        n_tot = normal_cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        t_tot = tumor_cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        tumor_logr.each_with_index do |l,i|
          if l[:count].to_i < 10 || normal_cov[i][:count].to_i < 10
            l[:_invalid] = true
            next
          end
          l[:name].gsub!(/[^\w]/,"")
          l[:logr] = Math.log((l[:count].to_f+10)/(normal_cov[i][:count].to_f+10)/(t_tot/n_tot),2).round(4)
        end
        tumor_logr.header[ tumor_logr.header.index(:count) ] = :logr

        tumor_logr.print config.tumor_exon_cnr

        normal_counts = exon_coverage normal_cov
        tumor_counts = exon_coverage tumor_cov

        cnr = {}

        File.open(config.tumor_gene_cnr,"w") do |f|
          f.puts "gene\tcounts"
          tumor_counts.each do |gene,counts|
            next if !normal_counts[gene] || normal_counts[gene] == 0
            cnr[gene] = Math.log(counts/normal_counts[gene],2).round(4)
            f.puts "#{gene}\t#{cnr[gene]}"
          end
        end

        muts = HashTable.new(config.insert_mutations)

        muts.header.insert( muts.header.index(:t_ref_count), :log_cnr )

        muts.each do |l|
          # fix bullshit
          gene = l[:"#gene"].split(/,/).first
          l[:log_cnr] = cnr[gene] || "-"
        end

        muts.print(config.tumor_mutations)
      end
    end
  end
end
