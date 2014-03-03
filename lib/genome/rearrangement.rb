#!/usr/bin/env ruby
require 'hash_table'
require 'vcf'
require 'fileutils'
module Genome
  class ExtractSclips
    include Pipeline::Step
    runs_on :tumor_samples, :chroms
    runs_tasks :extract_sclips

    class ExtractSclips
      include Pipeline::Task

      requires_file :sample_bam
      outs_file :sclip_file, :cover_file

      def run
        extract_sclips :bam => config.sample_bam, :chrom => config.chrom_name or error_exit "Could not run extractSclip!"
      end
    end
  end
  class CombineSclips
    include Pipeline::Step
    runs_on :tumor_samples
    runs_tasks :combine_sclips
    class CombineSclips
      include Pipeline::Task

      requires_file :chroms__cover_files
      outs_file :sample_cover_file

      def run
        run_cmd "cat #{config.chroms__cover_files.join(" ")} > #{config.sample_cover_file}" or error_exit "Could not combine coverage files!"
      end
    end
  end
  class RunCrest
    include Pipeline::Step
    runs_on :tumor_samples, :chroms
    runs_tasks :run_crest
    resources :threads => 12
    class RunCrest
      include Pipeline::Task

      requires_file :sample_bam, :normal_bam, :sample_cover_file
      outs_file :crest_rearr => :exists

      def run
        blat_check_status or error_exit "Couldn't connect to blat server"
        crest :tumor_bam => config.sample_bam, :normal_bam => config.normal_bam, :chrom => config.chrom_name, :sclip => config.sample_cover_file,
          :prefix => "#{config.sample_name}.#{config.chrom_name}" or error_exit "Could not run CREST!"
      end
    end
  end
  class CombineRearrs
    include Pipeline::Step
    runs_on :tumor_samples
    runs_tasks :combine_rearrs
    class CombineRearrs
      include Pipeline::Task

      requires_file :chroms__crest_rearrs => :exists
      outs_file :sample_crest_rearr

      def run
        run_cmd "cat #{config.chroms__crest_rearrs.join(" ")} > #{config.sample_crest_rearr}" or error_exit "Could not combine coverage files!"
      end
    end
  end
end
