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

  class CnvkitCoverage
    include Pipeline::Step
    runs_on :samples
    runs_tasks :target_cov, :anti_cov

    class TargetCov
      include Pipeline::Task

      requires_file :sample_bam
      outs_file :sample_target_cov
      def run
        cnvkit_coverage input: config.sample_bam, 
          target: config.cnvkit_target, output: config.sample_target_cov or error_exit "Could not run cnvkit"
      end
    end

    class AntiCov
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :sample_antitarget_cov
      def run
        cnvkit_coverage input: config.sample_bam,
          target: config.cnvkit_background, output: config.sample_antitarget_cov or error_exit "Could not run cnvkit"
      end
    end
  end

  class CnvkitMakeReference
    include Pipeline::Step
    runs_on :normal_samples
    runs_tasks :make_reference

    class MakeReference
      include Pipeline::Task
      requires_file :sample_target_cov, :sample_antitarget_cov
      outs_file :sample_reference_cnn
      def run
        cnvkit_reference inputs: [ config.sample_antitarget_cov, config.sample_target_cov ],
          output: config.sample_reference_cnn or error_exit "Could not run cnvkit"
      end
    end
  end

  class CnvkitFix
    include Pipeline::Step
    runs_on :tumor_samples
    runs_tasks :fix_cnr, :copy_seg

    class FixCnr
      include Pipeline::Task
      requires_file :sample_target_cov, :sample_antitarget_cov, :normal_reference_cnn
      outs_file :sample_cnr
      def run
        cnvkit_fix target: config.sample_target_cov, 
          background: config.sample_antitarget_cov, 
          reference: config.normal_reference_cnn, 
          output: config.sample_cnr or error_exit "Could not run cnvkit"
      end
    end

    class CopySeg
      include Pipeline::Task
      requires_file :sample_cnr
      outs_file :tumor_cnvkit_rdata, :tumor_cnvkit_seg

      def run
        # just pass these arguments to the R script
        r_script :segment, :doSegCbs, config.sample_cnr, config.tumor_cnvkit_rdata, config.tumor_cnvkit_seg, config.sample_name, config.segment_smoothing or error_exit "CBS segmentation failed"
      end
    end
  end
end
