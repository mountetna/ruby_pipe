# This script takes a bunch of samples (with replicates) and aligns them with tophat, then runs cufflinks on each of the samples, then runs cuffdiff to compare the samples to the parental
require 'pipeline'
require 'ribo/config'
require 'ribo/align'
require 'ribo/tophat'
require 'ribo/build'
require 'ribo/combine'
require 'ribo/coverage'
require 'ribo/qc'
require 'ribo/summary'
require 'ribo/align_rsem'
require 'ribo/babel'

module Ribo
  class ProfAlign 
    include Pipeline::Script
    runs_steps :align, :rsem_align, :tophat, :bwa_align, :combine, :coverage, :qc, :summary, :babel, :build

    def_module :default, :align => true,
      :tophat => true,
      :combine => true,
      :coverage => true,
      :qc => true,
      :summary => true

    def_module :setup, build: [ :create_transcript_model ]
    def_module :rsem, :align => [ :clip_fastq, :soak_ribo, :cull_non_ribo, :collect_rrna_metrics, :make_nonribo_fastq ],
      :rsem_align => true,
      :bwa_align => true,
      :coverage => true,
      :qc => true,
      :summary => true,
      :babel => true

    def_module :transcript_model_coverage, 
      coverage: [ :transcript_model_coverage ],
      summary: [ :summarize_transcript_model_cov, :summarize_qc ]

    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.first
        sample[:input_fastq] = args[1]
      end
      usage "input <sample name> <fastq>", "Input fastq file for the given sample"

      def qual_type(args)
        if args.first !~ /^(solexa|phred64)$/
          usage :qual_type
          return
        end
        config[:qual_type] = args.first
      end
      usage "qual_type <solexa|phred64>", "Specify the quality type for this set of fastq files"
    end
  end
end
