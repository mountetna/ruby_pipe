# This script takes a bunch of samples (with replicates) and aligns them with tophat, then runs cufflinks on each of the samples, then runs cuffdiff to compare the samples to the parental
require 'pipeline'
require 'rna/univ_geno'
require 'rna/tophat_align'
require 'rna/count_transcripts'
require 'rna/assemble_transcripts'
require 'rna/compare_expn'
require 'rna/univ_geno'
require 'rna/filter_muts'
require 'rna/qc'
require 'rna/config'

module Rna
  class PairedAlign 
    include Pipeline::Script
    runs_steps :tophat_align, :qc, :count_transcripts, :assemble_transcripts, :univ_geno, :filter_muts

    module Config
      include Pipeline::SampleConfig
      # this is config stuff that is particular to this sample

      def replicate_index
        case step
        when :tophat_align, :count_transcripts, :qc
          key_index_list(:replicates)[job_index]
        end
      end

      def sample_index
        case step
        # runs on every replicate
        when :tophat_align, :count_transcripts, :qc
          sample_index_list(:replicates)[job_index]
        # runs on every tumor sample
        when :filter_muts
          job_index+1
        end
      end

      def splits
        case step
        when :tophat_align, :count_transcripts, :qc
          replicates.size
        when :assemble_transcripts, :univ_geno
          1
        when :filter_muts
          samples.size - 1
        end
      end

      def scratch
        case step
        when :tophat_align, :count_transcripts, :qc
          replicate_scratch
        when :filter_muts
          sample_scratch
        when :univ_geno, :assemble_transcripts
          job_scratch
        end
      end
    end
  end
end
