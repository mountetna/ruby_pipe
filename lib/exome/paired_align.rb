require 'pipeline'
require 'exome/prep'
require 'exome/align'
require 'exome/dedup'
require 'exome/recal'
require 'exome/realign'
require 'exome/hybrid_qc'
require 'exome/mut_det'
require 'exome/mut_filter'
require 'exome/univ_geno'
require 'exome/config'
require 'exome/copy_number'
require 'exome/fastqc'
require 'exome/summarize'
require 'exome/snpeff'
require 'exome/somaticindel'

module Exome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :prep, 
      :fast_qc,
      :dump_fastqs, :combine_fastqs, :align, :dedup,
      :lane_recal, :table_recal,
      :patient_realign, 
      :make_samples,
      :hybrid_qc, :hybrid_qc_summary, 
      :sample_coverage, :copy_number, 
      :run_ascat,
      :mut_det, :mut_filter, :combine_muts,
      :univ_geno_call, :univ_geno_annotate,
      :run_absolute, :review_absolute,
      :summarize

    def_module :prep_pipe, {
      :prep => true
    }

    def_module :verify_fastqs, {
      :fast_qc => true
    }

    def_module :align_bams, {
      :dump_fastqs => true,
      :combine_fastqs => true,
      :align => true,
      :dedup => true
    }

    def_module :align_bwa_mem, {
      :align => [ :make_fastq_chunk, :align_mem, :verify_mate, :enforce_label ]
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
      :hybrid_qc => true,
      :hybrid_qc_summary => true,
    }

    def_module :ascat_purity, {
      :sample_coverage => true,
      :compute_normals => true,
      :copy_number => true,
      :run_ascat => true
    }

    def_module :absolute_purity, {
      :run_absolute => true,
      :review_absolute => true
    }

    def_module :compute_copy_number, {
      :prep => true,
      :sample_coverage => true,
      :copy_number => true
    }

    def_module :find_mutations, {
      :mut_det => true,
      :mut_filter => true,
      :combine_muts => true
    }

    def_module :find_mutations_indelocator, {
      :mut_det => [ :indelocator ]
    }

    def_module :find_mutations_pindel, {
      :mut_det => [ :pindel, :pindel_vcf, :patch_pindel_vcf ]
    }

    def_module :find_mutations_strelka, {
      :mut_det => [ :strelka ]
    }

    def_module :find_germline_muts, {
      :univ_geno_call => true,
      :univ_geno_annotate => true
    }

    def_module :find_mutations_somatic_indel_detector, {
      :mut_det => [ :mutect, :somatic_indel_detector, :patch_somatic_indel_vcf ],
      :mut_filter => [ :mutect_to_vcf, :snp_eff_annotate_mutect_vcf, :snp_eff_annotate_somatic_indel, :filter_muts_somatic_indel ],
      :combine_muts => true
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
      :find_mutations_somatic_indel_detector => true,
      :absolute_purity => true,
      :summarize => true
    }

    def exclude_task? task
      return nil
    end

    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.shift
        sample[:inputs] = []
        unless args.empty?
          sample[:inputs] += %x{ ls #{args.first} }.strip.split
        else
          getlines "your fastq files" do |l|
            sample[:inputs] += l.strip.split
          end
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
