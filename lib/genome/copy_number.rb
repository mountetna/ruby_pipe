#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
module Genome
  class CopyNumberPrep
    include Pipeline::Step
    runs_task :create_intervals_bed

    class CreateIntervalsBed
      include Pipeline::Task
      requires_file :interval_list
      dumps_file :interval_bed
      
      def run
        if !File.exists? config.interval_bed
          File.open(config.interval_bed,"w") do |f|
            File.foreach(config.interval_list) do |l|
              next if l =~ /^@/
              f.print l
            end
          end
        end
      end
    end
  end
  class CopyNumber
    include Pipeline::Step
    runs_tasks :compute_coverage, :compute_ratio, :copy_seg
    runs_on :tumor_samples

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
      requires_files :normal_cov, :tumor_cov
      outs_files :tumor_exon_cnr

      def exon_coverage(cov)
        tot = cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        cov.each_with_object(Hash.new(0)) do |l,m| 
          m[l[:name]] += (l[:count].to_i + 10) / tot 
        end
      end

      def run
        header = [ :chr, :start, :stop, :strand, :name, :count ]
        normal_cov = HashTable.new(config.normal_cov, :header => header)
        tumor_cov = HashTable.new(config.tumor_cov, :header => header)
        tumor_logr = HashTable.new(config.tumor_cov, :header => header)

        n_tot = normal_cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        t_tot = tumor_cov.inject(0) { |m,l| m += l[:count].to_i + 10 }.to_f
        tumor_logr.each_with_index do |l,i|
          if l[:count].to_i < 10 || normal_cov[i][:count].to_i < 10
            l[:_invalid] = true
            next
          end
          l[:logr] = Math.log((l[:count].to_f+10)/(normal_cov[i][:count].to_f+10)/(t_tot/n_tot),2).round(4)
          l[:normal_count] = normal_cov[i][:count]
          l[:tumor_count] = l[:count]
        end
        tumor_logr.header[ tumor_logr.header.index(:count) ] = [ :tumor_count, :normal_count, :logr ]
        tumor_logr.header.flatten!

        tumor_logr.print config.tumor_exon_cnr
      end
    end

    class CopySeg
      include Pipeline::Task
      requires_file :tumor_exon_cnr
      outs_file :tumor_cnr_rdata, :tumor_cnr_seg

      def run
        # just pass these arguments to the R script
        r_script :segment, config.tumor_exon_cnr, config.tumor_cnr_rdata, config.tumor_cnr_seg, config.sample_name
      end
    end
    class ComputePurityPloidy
      include Pipeline::Task
      requires_file :tumor_cnr_seg

      def run
        r_script :absolute, config.sample_name, config.tumor_cnr_seg, config.tumor_maf, config.absolute_scratch
      end
    end
  end
end
