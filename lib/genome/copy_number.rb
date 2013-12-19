#!/usr/bin/env ruby
require 'hash_table'
require 'vcf'
require 'fileutils'
module Genome
  class CopyNumber
    include Pipeline::Step
    runs_tasks :compute_coverage, :compute_ratio, :copy_seg_cbs, :compute_purity_ploidy
    has_tasks :compute_baf, :compute_coverage, :compute_ratio, :copy_seg_cbs, :copy_seg_pscbs, :compute_purity_ploidy
    runs_on :tumor_samples
    resources :threads => 1, :walltime => 50

    class ComputeBaf
      include Pipeline::Task
      requires_file :ug_filtered_vcf 
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

      def run
        header = [ :chr, :start, :stop, :count ]
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

    class CopySegCbs
      include Pipeline::Task
      requires_file :tumor_exon_cnr
      outs_file :tumor_cnr_rdata, :tumor_cnr_seg

      def run
        r_script :segment, :doSegCbs, config.tumor_exon_cnr, config.tumor_cnr_rdata, config.tumor_cnr_seg, config.sample_name
      end
    end
    class CopySegPscbs
      include Pipeline::Task
      requires_file :tumor_ratio, :normal_baf, :tumor_baf
      outs_file :tumor_cnr_rdata_pscbs, :tumor_cnr_seg_pscbs

      def run
        # just pass these arguments to the R script
        r_script :segment, :doSegPscbs, config.tumor_ratio, config.tumor_baf, config.normal_baf, config.tumor_cnr_rdata_pscbs or error_exit "CBS segmentation failed"
       end
    end
    class ComputePurityPloidy
      include Pipeline::Task
      requires_file :tumor_cnr_seg, :all_muts_maf
      outs_file :absolute_rdata

      def run
        r_script :absolute, :callSample, config.sample_name, config.tumor_cnr_seg, config.all_muts_maf, config.absolute_scratch or error_exit "Absolute failed"
      end
    end
  end

  class Absolute 
    include Pipeline::Step
    runs_tasks :compute_purity_ploidy
    runs_on :tumor_samples

    class ComputePurityPloidy
      include Pipeline::Task
      requires_file :tumor_cnr_seg, :all_muts_maf
      outs_file :absolute_rdata

      def run
        r_script :absolute, :callSample, config.sample_name, config.tumor_cnr_seg, config.all_muts_maf, config.absolute_scratch or error_exit "Absolute failed"
      end
    end
  end
 
  class ReviewAbsolute
    include Pipeline::Step
    runs_tasks :create_review_object, :extract_review_results

    class CreateReviewObject
      include Pipeline::Task

      requires_files :absolute_rdatas
      dumps_file :review_table

      def run
        r_script :absolute, :createReview, config.cohort_name, config.absolute_review_dir, *config.absolute_rdatas
      end
    end
    class ExtractReviewResults
      include Pipeline::Task

      requires_file :reviewed_table
      def run
        r_script :absolute, :extractReview, config.reviewed_table, "Genome.Pipeline", config.absolute_modes, config.absolute_review_dir, config.cohort_name
      end
    end
  end
end
