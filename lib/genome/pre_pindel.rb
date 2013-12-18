#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'vcf'

module Genome
  class PrePindel 
    include Pipeline::Step
    runs_tasks :make_pre_pindel 
    resources :threads => 1, :walltime => 50
    runs_on :samples

    class MakePrePindel 
      include Pipeline::Task
      requires_files :sample_bam 
      outs_file :temp_output_for_pindel

      def run
        log_info "Creating files for Pindel"
        run_cmd "#{config.samtools_dir}/samtools view #{config.sample_bam} | sam2pindel - #{config.temp_output_for_pindel} 300 #{config.sample_name} 0 Illumina-PairEnd" 
      end
    end 
  end

  class PrePindelTwo
    include Pipeline::Step
    runs_tasks :pre_pindel_cat
    resources :threads => 1, :walltime => 50
    runs_on :tumor_samples
    
    class PrePindelCat
      include Pipeline::Task
      requires_files :temp_output_for_pindels
      outs_file :output_for_pindel

      def run
        log_info "Merging files necessary for Pindel"
        run_cmd "cat #{config.temp_output_for_pindel(config.normal)} #{config.temp_output_for_pindel(config.sample)} > #{config.output_for_pindel}"
      end
    end
  end

end
