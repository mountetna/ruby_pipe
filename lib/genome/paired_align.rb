require 'pipeline'
require 'genome/fix_mate'
require 'genome/recal'
require 'genome/hybrid_qc'
require 'genome/mut_det'
require 'genome/config'
require 'genome/copy_number'

module Genome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :fix_mate, :library_merge, :recal, :library_split, :make_samples, :hybrid_qc, :hybrid_qc_summary, :copy_number, :mut_det, :mut_filter

    def_module :create_bams, {
      :fix_mate => true,
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
      :copy_number => true,
    }

    def_module :find_mutations, {
      :mut_det => true,
      :mut_filter => true
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