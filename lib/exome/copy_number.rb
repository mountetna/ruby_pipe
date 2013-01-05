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
      outs_files :tumor_cnr, :tumor_mutations

      def exon_coverage(cov)
        tot = cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        cov.each_with_object(Hash.new(0)) do |l,m| 
          m[l[:name]] += (l[:count].to_i + 10) / tot 
        end
      end

      def run
        normal_cov = HashTable.new(config.normal_cov, [ :chr, :start, :stop, :strand, :name, :count ])
        tumor_cov = HashTable.new(config.tumor_cov, [ :chr, :start, :stop, :strand, :name, :count ])

        normal_counts = exon_coverage normal_cov
        tumor_counts = exon_coverage tumor_cov

        cnr = {}

        File.open(config.tumor_cnr,"w") do |f|
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
