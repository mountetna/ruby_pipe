#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

module Exome
  class Summarize
    include Pipeline::Step
    runs_tasks :concat_mafs, :concat_segs, :combine_absolute_pdfs, :copy_main_log, :copy_metrics_summary
    runs_on :cohort

    class ConcatMafs
      include Pipeline::Task
      requires_files :tumor_somatic_mafs
      outs_files :combined_somatic_maf

      def run
        @somatic_maf = Maf.new

        config.tumor_samples.each do |sample|
          m = Maf.read config.tumor_maf(sample)
          @somatic_maf.headers = m.headers
          @somatic_maf.lines.concat m.lines
        end
        @somatic_maf.write config.combined_somatic_maf
      end
    end
    class ConcatSegs
      include Pipeline::Task
      requires_files :tumor_cnr_segs
      outs_files :combined_seg

      def run
        @seg = HashTable.new nil

        config.tumor_samples.each do |sample|
          sseg = HashTable.new config.tumor_cnr_seg(sample)
          @seg.header = sseg.header
          sseg.each do |line|
            @seg.add_line line
          end
        end
        @seg.write config.combined_seg
      end
    end
    class CombineAbsolutePdfs
      include Pipeline::Task
      requires_files :tumor_absolute_pdfs
      outs_files :combined_absolute_pdf

      def run
        merge_pdfs config.tumor_samples.map{|s| config.absolute_pdf(s) }, config.combined_absolute_pdf or error_exit "Merge_pdfs failed"
      end
    end
    class CopyMainLog
      include Pipeline::Task
      requires_file :main_log
      outs_file :main_log_copy

      def run
        FileUtils.cp config.main_log, config.main_log_copy
      end
    end
    class CopyMetricsSummary
      include Pipeline::Task
      requires_file :qc_summary
      outs_file :qc_summary_copy

      def run
        FileUtils.cp config.qc_summary, config.qc_summary_copy
      end
    end
  end
end
