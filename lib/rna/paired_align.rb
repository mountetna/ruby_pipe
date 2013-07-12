# This script takes a bunch of samples (with replicates) and aligns them with tophat, then runs cufflinks on each of the samples, then runs cuffdiff to compare the samples to the parental
require 'pipeline'
require 'rna/univ_geno'
require 'rna/tophat_align'
require 'rna/count_transcripts'
require 'rna/assemble_transcripts'
require 'rna/compare_expn'
require 'rna/univ_geno'
require 'rna/diff_exp'
require 'rna/filter_muts'
require 'rna/qc'
require 'rna/config'

module Rna
  class PairedAlign 
    include Pipeline::Script
    runs_steps :tophat_align, :qc, :count_transcripts, :diff_exp, :assemble_transcripts #, :univ_geno, :filter_muts

    def_module :rsem, :tophat_align => true,
      :qc => true,
      :count_transcripts => [ :rsem_count ]

    def_module :default, :tophat_align => true,
      :qc => true,
      :count_transcripts => true,
      :assemble_transcripts => true,
      :diff_exp => true

    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.first
        inputs = []
        getlines "your fastq files" do |l|
          inputs += l.strip.split
        end
        sample[:replicates] ||= []
        sample[:replicates].push :inputs => pair_fastqs(inputs)
      end
      usage "input <sample name>", "Input fastqs for the given sample"

      def diff_exp(args)
        sample = find_sample args[0]
        normal = find_sample args[1]
        if !sample || !normal
          puts "No sample of that name.".red
          return
        end
        if sample == normal
          puts "Sample and normal are the same".red
          return
        end

        sample[:diff_exp] ||= []
        sample[:diff_exp].push :normal_name => normal[:sample_name]
      end
      usage "diff_exp <sample name> <normal sample name>", "Differential expression between two samples"

      def frag_size(args)
        config[:frag_size] = args.first
      end
      usage "frag_size <size in bp>", "Set the fragment size for the given sample (reads+insert)."
    end
  end
end
