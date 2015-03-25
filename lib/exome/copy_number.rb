#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'vcf'

module Exome

  class SampleCoverage
    include Pipeline::Step
    runs_tasks :compute_coverage # :correct_gc
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
    
    class CorrectGc
      include Pipeline::Task
      requires_files :sample_cov
      dumps_files :sample_cov_gc

      def run
        r_script :segment, :doGcCorrect, config.sample_cov, config.sample_cov_gc, config.reference_gc or error_exit "GC correction failed for tumor!"
      end
    end
  end

  class CoverageTable < HashTable
    print_columns
    columns :seqname, :start, :stop, :strand, :name, :count
    types start: :int, stop: :int, count: :int
  end

  class CopyRatioTable < HashTable
    print_columns
    columns :chrom, :start, :stop, :strand, :name, :normal_count, :tumor_count, :logr

    class CopyRatio < HashTable::Row
      def logr
         Math.log(pseudo(tumor_count)/pseudo(normal_count) / (@table.tumor_reads/@table.normal_reads),2).round(4)
      end

      def enough_reads?
        normal_count > 10
      end

      def pseudo(count)
        count.to_f + 10
      end
    end

    def initialize normal_cov, sample_cov
      super []
      sample_cov.lazy.zip(normal_cov) do |scov,ncov|
        self << {
          chrom: ncov.seqname,
          start: ncov.start,
          stop: ncov.stop,
          strand: ncov.strand,
          name: ncov.name,
          normal_count: ncov.count,
          tumor_count: scov.count
        }
      end
    end

    def tumor_reads
      @tumor_reads ||= total_reads :tumor_count
    end

    def normal_reads
      @normal_reads ||= total_reads :normal_count
    end

    def total_reads type
      inject(0) do |sum,line|
        sum += line.pseudo(line.send(type)) if line.enough_reads?
        sum
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
        normal_cov = CoverageTable.new
        normal_cov.parse config.normal_cov
        sample_cov = CoverageTable.new
        sample_cov.parse config.sample_cov
        tumor_logr = CopyRatioTable.new normal_cov, sample_cov

        tumor_logr.print config.sample_exon_cnr do |l|
          l.enough_reads? ? l : nil
        end
      end
    end

    class CopySeg
      include Pipeline::Task
      requires_file :sample_exon_cnr
      outs_file :tumor_cnr_rdata, :tumor_cnr_seg

      def run
        # just pass these arguments to the R script
        r_script :segment, :doSegCbs, config.sample_exon_cnr, config.tumor_cnr_rdata, config.tumor_cnr_seg, config.sample_name, config.segment_smoothing or error_exit "CBS segmentation failed"
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
        tb = HashTable.new :columns => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        nb = HashTable.new :columns => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        vcf.each do |v|
          next if !v.genotype(config.normal.sample_id).heterozygous?
          next if !v.genotype(config.sample.sample_id).callable?
          next if v.filter != "PASS"
          next if v.id[0..1] != "rs"
          tb.add_line :chromosome => v.seqname, :position => v.pos,
            :BAF => v.genotype(config.sample.sample_id).alt_freq,
            :Alt_count => v.genotype(config.sample.sample_id).alt_count,
            :Tot_count => v.genotype(config.sample.sample_id).depth
          nb.add_line :chromosome => v.seqname, :position => v.pos,
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
      outs_file :absolute_rdata, :absolute_pdf

      def run
        r_script :absolute, :callSample, config.sample_name, config.tumor_cnr_seg, config.all_muts_maf, config.absolute_scratch or error_exit "Absolute failed"

        FileUtils.cp config.absolute_scratch_pdf, config.absolute_pdf if File.exists? config.absolute_scratch_pdf
        FileUtils.cp config.absolute_scratch_rdata, config.absolute_rdata if File.exists? config.absolute_scratch_rdata
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
        r_script :absolute, :createReview, config.cohort_name, config.absolute_review_dir, config.absolute_rdatas.join(",") or error_exit "Absolute failed"
      end
    end
    class ExtractReviewResults
      include Pipeline::Task

      requires_file :reviewed_table, :absolute_modes
      outs_file  :absolute_calls

      def run
        r_script :absolute, :extractReview, config.reviewed_table, "Exome.Pipeline", config.absolute_modes, config.absolute_review_dir, config.cohort_name or error_exit "Absolute failed"
        FileUtils.mv config.absolute_calls_scratch, config.absolute_calls 
        config.tumor_samples.each do |sample|
          FileUtils.mv config.absolute_segs_scratch(sample), config.tumor_absolute_seg(sample)
        end
      end
    end
  end
end
