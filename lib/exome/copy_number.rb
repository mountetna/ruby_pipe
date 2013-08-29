#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'vcf'

module Exome
  class CoverageTable < HashTable
    header_on
    class Coverage < HashLine
      def p_count
        count.to_f + 10
      end
      def low_count?
        count.to_f < 10
      end
    end
    line_class :coverage


    def total
      @total ||= inject(0) do |sum,line|
        sum += line.p_count
      end.to_f
    end

    def initialize file
      super file, :header => [ :chr, :start, :stop, :strand, :name, :count ]
    end
  end

  class SampleCoverage
    include Pipeline::Step
    runs_tasks :compute_coverage
    runs_on :samples
    audit_report :sample_name

    class ComputeCoverage
      include Pipeline::Task
      requires_files :sample_bam, :interval_bed
      dumps_files :sample_cov

      def run
        coverage_bed config.sample_bam, config.interval_bed, config.sample_cov or error_exit "Computing coverage failed."
      end
    end
  end

  class PrepNormal
    include Pipeline::Step
    runs_task :compute_total_coverage
    
    class ComputeTotalCoverage
      include Pipeline::Task
      requires_files :normal_sample_covs
      dumps_file :total_normal_cov

      def run
        total_cov = CoverageTable.new config.normal_sample_covs.first
        normal_covs = config.normal_sample_covs.map{|f| CoverageTable.new f }
        total_cov.each_with_index do |line,i|
          line[:count] = normal_covs.inject(0) do |sum,nc| 
            sum += nc[i].count.to_f
          end
        end
        total_cov.print config.total_normal_cov
      end
    end
  end

  class ComputeNormals
    include Pipeline::Step
    runs_task :compute_normal_logr, :copy_seg
    runs_on :normal_samples

    class ComputeNormalLogr
      include Pipeline::Task
      requires_file :sample_cov, :total_normal_cov
      outs_files :sample_exon_cnr

      def run
        total_cov = CoverageTable.new config.total_normal_cov
        normal_logr = CoverageTable.new config.sample_cov
        normal_total = normal_logr.total

        normal_logr.each_with_index do |line,i|
          line[:total_count] = total_cov[i].count
          line[:normal_count] = line.count
          if line.low_count?
            line[:logr] = "NA"
          else
            line[:logr] = Math.log(line.p_count/total_cov[i].p_count/(normal_total/total_cov.total),2).round(4)
          end
        end
        normal_logr.header[ normal_logr.header.index(:count) ] = [ :normal_count, :total_count, :logr ]
        normal_logr.header.flatten!

        normal_logr.print config.sample_exon_cnr
      end
    end

    class CopySeg
      include Pipeline::Task
      requires_file :sample_exon_cnr
      outs_file :tumor_cnr_rdata, :tumor_cnr_seg

      def run
        # just pass these arguments to the R script
        r_script :segment, :doSegCbs, config.sample_exon_cnr, config.tumor_cnr_rdata, config.tumor_cnr_seg, config.sample_name or error_exit "CBS segmentation failed"
      end
    end
  end

  class CopyNumber
    include Pipeline::Step
    runs_tasks :compute_ratio, :copy_seg #:absolute_purity_ploidy
    runs_on :tumor_samples
    audit_report :sample_name, :normal_name

    class ComputeRatio
      include Pipeline::Task
      requires_files :normal_cov, :sample_cov
      outs_files :sample_exon_cnr

      def run
        normal_cov = CoverageTable.new config.normal_cov
        tumor_logr = CoverageTable.new config.sample_cov
        total = tumor_logr.total

        tumor_logr.each_with_index do |line,i|
          line[:normal_count] = normal_cov[i].count
          line[:tumor_count] = line.count
          if line.low_count? || normal_cov[i].low_count?
            line[:logr] = "NA"
          else
            line[:logr] = Math.log(line.p_count/normal_cov[i].p_count/(total/normal_cov.total),2).round(4)
          end
        end
        tumor_logr.header[ tumor_logr.header.index(:count) ] = [ :tumor_count, :normal_count, :logr ]
        tumor_logr.header.flatten!

        tumor_logr.print config.sample_exon_cnr
      end
    end

    class CopySeg
      include Pipeline::Task
      requires_file :sample_exon_cnr
      outs_file :tumor_cnr_rdata, :tumor_cnr_seg

      def run
        # just pass these arguments to the R script
        r_script :segment, :doSegCbs, config.sample_exon_cnr, config.tumor_cnr_rdata, config.tumor_cnr_seg, config.sample_name or error_exit "CBS segmentation failed"
      end
    end

  end

  class RunAscat
    include Pipeline::Step
    runs_on :tumor_samples
    audit_report :sample_name, :normal_name
    runs_tasks :compute_baf, :ascat_purity_ploidy

    class ComputeBaf
      include Pipeline::Task
      requires_file :ug_snps_vcf
      outs_file :tumor_baf, :normal_baf

      def run
        vcf = VCF.read config.ug_snps_vcf
        tb = HashTable.new nil, :header => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        nb = HashTable.new nil, :header => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        vcf.each do |v|
          next if !v.genotype(config.normal.sample_id).heterozygous?
          next if !v.genotype(config.sample.sample_id).callable?
          next if v.filter != "PASS"
          next if v.id[0..1] != "rs"
          tb.add_line :chromosome => v.chrom, :position => v.pos,
            :BAF => v.genotype(config.sample.sample_id).alt_freq,
            :Alt_count => v.genotype(config.sample.sample_id).alt_count,
            :Tot_count => v.genotype(config.sample.sample_id).depth
          nb.add_line :chromosome => v.chrom, :position => v.pos,
            :BAF => v.genotype(config.normal.sample_id).alt_freq,
            :Alt_count => v.genotype(config.normal.sample_id).alt_count,
            :Tot_count => v.genotype(config.normal.sample_id).depth
        end
        tb.print config.tumor_baf
        nb.print config.normal_baf
      end
    end

    class AscatPurityPloidy
      include Pipeline::Task
      requires_file :sample_exon_cnr, :tumor_baf, :normal_exon_cnr, :normal_baf, :interval_bed
      outs_file :tumor_ascat_rdata, :tumor_ascat_txt

      def run
        r_script :segment, :doAscatPurityPloidy, config.tumor_baf, config.normal_baf, config.sample_exon_cnr, config.normal_exon_cnr, config.interval_bed, config.tumor_ascat_rdata, config.tumor_ascat_txt or error_exit "ASCAT failed"
      end
    end
  end

  class RunAbsolute
    include Pipeline::Step
    runs_on :tumor_samples
    audit_report :sample_name, :normal_name
    runs_tasks :absolute_purity_ploidy
    class AbsolutePurityPloidy
      include Pipeline::Task
      requires_file :sample_exon_cnr, :all_muts_maf
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
        r_script :absolute, :createReview, config.cohort_name, config.absolute_review_dir, *config.absolute_rdatas or error_exit "Absolute failed"
      end
    end
    class ExtractReviewResults
      include Pipeline::Task

      requires_file :reviewed_table
      def run
        r_script :absolute, :extractReview, config.reviewed_table, "Exome.Pipeline", config.absolute_modes, config.absolute_review_dir, config.cohort_name or error_exit "Absolute failed"
      end
    end
  end
end
