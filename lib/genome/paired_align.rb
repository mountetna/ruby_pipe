require 'pipeline'
require 'genome/align'
require 'genome/sam_merge'
require 'genome/fix_mate'
require 'genome/recal'
require 'genome/collect_qc'
require 'genome/mut_det'
require 'genome/config'
require 'genome/copy_number'
require 'genome/indel'

module Genome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :dump_fastqs, :combine_fastqs, :align,
      :merge,
      :lane_recal, :table_recal,
      :patient_realign, 
      :make_samples, 
      :collect_qc, :collect_qc_summary,
      :mut_det, :indel_det, :copy_number, 
      :variant_det, :mut_filter, :review_absolute, :merge_snp

    def_module :align_bams, {
      :dump_fastqs => true,
      :combine_fastqs => true,
      :align => true,
      :merge => true
    }

    def_module :recal_by_lane, {
      :lane_recal => true,
      :table_recal => true
    }

    def_module :realign_by_patient, {
      :patient_realign => true,
      :patient_split => true
    }

    def_module :finalize_samples, {
      :make_samples => true
    }

    def_module :create_bams, {
      :align_bams => true,
      :recal_by_lane => true,
      :realign_by_patient => true,
      :finalize_samples => true
    }

    def_module :calculate_qc, {
      :collect_qc => true,
      :collect_qc_summary => true,
    }

    def_module :compute_copy_number, {
      :copy_number => true,
      :merge_snp => true,
      :format_rdata => true,
      :absolute => true,
      :review_absolute => true,
    }

    def_module :find_mutations, {
      :variant_det => true,
      :merge_snp => true,
      :mut_det => true,
      #:indel_det => true,
      :mut_filter => true,
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
