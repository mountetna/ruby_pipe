require 'pipeline'
require 'exome/align'
require 'exome/merge'
require 'exome/recal'
require 'exome/hybrid_qc'
require 'exome/mut_det'
require 'exome/config'
require 'exome/copy_number'

module Exome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :align, :merge, :library_merge, :recal, :library_split, :make_samples, :hybrid_qc, :hybrid_qc_summary, :copy_number_prep, :copy_number, :mut_det, :mut_filter

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
