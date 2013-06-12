require 'pipeline'
require 'exome/align'
require 'exome/merge'
require 'exome/recal'
require 'exome/hybrid_qc'
require 'exome/mut_det'
require 'exome/univ_geno'
require 'exome/config'
require 'exome/copy_number'

module Exome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :align, :merge, :library_merge, :recal, :library_split, :make_samples, :hybrid_qc, :hybrid_qc_summary, :copy_number_prep, :copy_number, :review_absolute, :mut_det, :mut_filter, :univ_geno_normals

    def_module :create_bams, {
      :align => true,
      :merge => true,
      :library_merge => true,
      :recal => true,
      :library_split => true,
      :make_samples => true, 
    }

    def_module :calculate_qc, {
      :hybrid_qc => true,
      :hybrid_qc_summary => true,
    }

    def_module :compute_copy_number, {
      :copy_number_prep => true,
      :copy_number => true,
      :review_absolute => true
    }

    def_module :find_mutations, {
      :mut_det => true,
      :mut_filter => true
    }

    def_module :find_normal_mutations, {
      :univ_geno_normals => true
    }

    def_module :mut_filter_annovar, {
      :mut_filter => [ :concat_chroms, :filter_muts_annovar ]
    }

    def_module :default, {
      # the default sequence of events. The order is dictated by runs_steps
      :create_bams => true,
      :calculate_qc => true,
      :compute_copy_number => true,
      :find_mutations => true
    }

    def exclude_task? task
      return nil
    end

    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.first
        sample[:inputs] = []
        getlines "your fastq files" do |l|
          sample[:inputs] += l.strip.split
        end
        sample[:inputs] = pair_fastqs sample[:inputs]
      end
      usage "input <sample name>", "Input fastqs for the given sample"

      def intervals(args)
        config[:interval_list] = args.first
      end
      usage "intervals <interval list>", "Set the interval list file for the given sample."

      def frag_size(args)
        config[:frag_size] = args.first
      end
      usage "frag_size <size in bp>", "Set the fragment size for the given sample (reads+insert)."

      def normal(args)
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

        sample[:normal_name] = normal[:sample_name]
      end
      usage "normal <sample name> <normal sample name>", "Set the normal name for the given sample (the first sample is assumed to be the normal)"
    end
  end
end
