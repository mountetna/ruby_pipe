require 'pipeline'
require 'rna/univ_geno'
require 'rna/tophat_align'
require 'rna/count_transcripts'
require 'rna/assemble_transcripts'
require 'rna/compare_expn'
require 'rna/univ_geno'
require 'rna/diff_exp'
require 'rna/splice_count'
require 'rna/detect_fusions'
require 'rna/filter_muts'
require 'rna/qc'
require 'rna/config'

module Rna
  class PairedAlign 
    include Pipeline::Script
    runs_steps :rsem_count, :rsem_single_count, :rsem_format, :tophat_align, :qc, :qc_summary, :cufflinks_count, :cuff_diff_exp, :assemble_transcripts, :assemble_rsem_transcripts, :deseq_diff_exp, :splice_count, :detect_fusions #, :univ_geno, :filter_muts

    def_module :rsem,
      :rsem_count => true,
      :rsem_format => true,
      :qc => true,
      :qc_summary => true,
      :assemble_rsem_transcripts => true,
      :deseq_diff_exp => true

    def_module :rsem_single_end,
      rsem_count: [ :rsem_single_count ]

    def_module :count_splice, :splice_count => true

    def_module :fusion_detect, :detect_fusions => true

    def_module :cufflinks_count_denovo, :cufflinks_count => [ :cufflink_denovo, :format_transcript, :sort_sam, :count_coverage ]

    def_module :default, :tophat_align => true,
      :qc => true,
      :cufflinks_count => true,
      :assemble_transcripts => true,
      :cuff_diff_exp => true

    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.shift
        inputs = []
        unless args.empty?
          inputs += %x{ ls #{args.first} }.strip.split
        else
          getlines "your fastq files" do |l|
            inputs += l.strip.split
          end
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
