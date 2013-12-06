#!/usr/bin/env ruby
require '/taylorlab/lib/ruby/hash_table'
require 'fileutils'
require '/home/changmt/lib/ruby/vcf'

module Genome
  class CopyNumber2
    include Pipeline::Step
    runs_tasks :compute_baf, :copy_seg
    resources :threads => 1, :walltime => 50
    runs_on :tumor_samples
=begin
    class ComputeCoverage
      include Pipeline::Task
      requires_files :tumor_bam, :normal_bam, :interval_bed
      dumps_files :tumor_cov, :normal_cov

      def run
        coverage_bed config.normal_bam, config.interval_bed, config.normal_cov or error_exit "Computing normal coverage failed."
        coverage_bed config.tumor_bam, config.interval_bed, config.tumor_cov or error_exit "Computing tumor coverage failed."
      end
    end
=end
    class ComputeBaf
      include Pipeline::Task
      requires_file #:ug_filtered_vcf 
      outs_file :tumor_baf, :normal_baf

      def run
        vcf = VCF.read config.ug_filtered_vcf
        tb = HashTable.new nil, :header => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        nb = HashTable.new nil, :header => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        vcf.each do |v|
          next if !v.genotype(config.normal_name).heterozygous?
          next if !v.genotype(config.sample_name).callable?
          next if v.filter != "PASS"
          next if v.id[0..1] != "rs"
          tb.add_line :chromosome => v.chrom, :position => v.pos,
            :BAF => v.genotype(config.sample_name).alt_freq,
            :Alt_count => v.genotype(config.sample_name).alt_count,
            :Tot_count => v.genotype(config.sample_name).depth
          nb.add_line :chromosome => v.chrom, :position => v.pos,
            :BAF => v.genotype(config.normal_name).alt_freq,
            :Alt_count => v.genotype(config.normal_name).alt_count,
            :Tot_count => v.genotype(config.normal_name).depth
        end
        tb.print config.tumor_baf
        nb.print config.normal_baf
      end
    end    

    class CopySeg
      include Pipeline::Task
      requires_file :tumor_ratio, :normal_baf, :tumor_baf
      outs_file :tumor_cnr_rdata_pscbs, :tumor_cnr_seg_pscbs

      def run
        # just pass these arguments to the R script
        r_script :segment, :doSegPscbs, config.tumor_ratio, config.tumor_baf, config.normal_baf, config.tumor_cnr_rdata_pscbs or error_exit "CBS segmentation failed"
       end
    end

  end
end
