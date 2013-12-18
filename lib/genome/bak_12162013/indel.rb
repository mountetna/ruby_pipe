#!/usr/bin/env ruby

require '/home/changmt/lib/ruby/hash_table'
require 'fileutils'
require '/home/changmt/lib/ruby/mutect'
require '/home/changmt/lib/ruby/vcf'
require '/home/changmt/lib/ruby/maf'

module Genome
  class IndelDet
    include Pipeline::Step
    #runs_tasks :somatic_indels, :annotate_indels
    runs_tasks :indel_locator
    runs_on :tumor_samples, :chroms

    class IndelLocator
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      outs_files :indelocator_bed, :indelocator_metrics

      def run
        indelocator :somatic => true, :"input_file:normal" => config.normal_bam, :"input_file:tumor" => config.tumor_bam, :verboseOutput => config.indelocator_output, :metrics_file => config.indelocator_metrics, :out => config.indelocator_bed, :intervals => config.chrom or error_exit "Indelocator failed"
      end
    end
  end
=begin
    class SomaticIndels 
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :interval_bed
      outs_file :raw_indel_vcf

      def run
        log_info "Running Somatic Indel Detector..."
        gatk :somatic_indel_detector, "input_file:normal" => config.normal_bam,
                "input_file:tumor" => config.tumor_bam,
                "intervals" => config.interval_bed,
                "maxNumberOfReads" => 10000,
                "window_size" => 225,
                "filter_expressions" => '"N_COV<8||T_COV<14||T_INDEL_F<0.1||T_INDEL_CF<0.7"',
                "out" => config.raw_indel_vcf or error_exit "Indel detection failed"
      end
    end
    class AnnotateIndels
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :raw_indel_vcf, :interval_bed
      outs_file :annot_indel_vcf

      def run
        log_info "Annotating raw indel calls..."
        gatk :variant_annotator,
                :variant => config.raw_indel_vcf,
                :intervals => config.interval_bed,
                :"input_file:normal" => config.normal_bam,
                :"input_file:tumor" => config.tumor_bam,
                :dbsnp => config.dbsnp_vcf,
                :group => "StandardAnnotation",
                :out => config.annot_indel_vcf or error_exit "Indel annotation failed"
      end
    end
  end
=end
end
