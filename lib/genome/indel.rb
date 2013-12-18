#!/usr/bin/env ruby

require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

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
end
