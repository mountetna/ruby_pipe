# This script takes a bunch of samples (with replicates) and aligns them with tophat, then runs cufflinks on each of the samples, then runs cuffdiff to compare the samples to the parental
require 'pipeline'
require 'ribo/config'
require 'ribo/align'
require 'ribo/tophat'
require 'ribo/combine'
require 'ribo/coverage'
require 'ribo/qc'

module Ribo
  class ProfAlign 
    include Pipeline::Script
    runs_steps :align, :tophat, :combine, :coverage, :qc
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
