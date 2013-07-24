require 'pipeline'
require 'exome/prep'
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
    runs_steps :prep, 
      :dump_fastqs, :combine_fastqs, :align, :merge,
      :lane_merge, :lane_recal, :lane_table_recal, :lane_split,
      :patient_merge, :patient_realign, :patient_split, 
      :make_samples,
      :hybrid_qc, :hybrid_qc_summary, :copy_number, 
      :mut_det, :univ_geno_normals, :mut_filter, :review_absolute

    def_module :prep_pipe, {
      :prep => true
    }

    def_module :align_bams, {
      :dump_fastqs => true,
      :combine_fastqs => true,
      :align => true,
      :merge => true
    }

    def_module :recal_by_lane, {
      :lane_merge => true,
      :lane_recal => true,
      :lane_table_recal => true,
      :lane_split => true
    }

    def_module :realign_by_patient, {
      :patient_merge => true,
      :patient_realign => true,
      :patient_split => true
    }

    def_module :create_bams, {
      :align_bams => true,
      :recal_by_lane => true,
      :realign_by_patient => true,
      :make_samples => true
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

    def_module :find_mutations_indelocator, {
      :mut_det => [ :indelocator ]
    }

    def_module :find_normal_mutations, {
      :univ_geno_normals => true
    }

    def_module :mut_filter_annovar, {
      :mut_filter => [ :concat_chroms, :filter_muts_annovar ]
    }

    def_module :default, {
      # the default sequence of events. The order is dictated by runs_steps
      :prep_pipe => true,
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
