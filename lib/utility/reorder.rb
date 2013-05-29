require 'pipeline'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def_var :input_bam do |s| (s || sample).input_bam end

    dir_tree({
      ":output_dir" => {
        "@sample_name.reordered.bam" => :output_bam
      }
    })
  end

  class ReorderBam
    include Pipeline::Step
    runs_tasks :reorder_bam
    runs_on :samples

    class ReorderBam
      include Pipeline::Task
      requires_file :input_bam
      outs_file :output_bam

      def run
        log_info "Reordering bam file according to dictionary"
        picard :reorder_sam, :I => config.input_bam, :O => config.output_bam, :ALLOW_INCOMPLETE_DICT_CONCORDANCE => :true, :REFERENCE => config.reference_fa or error_exit "Could not reorder bam."
      end
    end
  end
  class Reorder
    include Pipeline::Script
    runs_steps :reorder_bam
    def_module :default, :reorder_bam => true
    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input_bam(args)
        sample = get_sample args.first
        sample[:input_bam] = args[1]
      end
      usage "input_bam <sample name> <input_bam>", "Input bam file for the given sample"

      def output_bam(args)
        sample = get_sample args.first
        sample[:output_bam] = args[1]
      end
      usage "output_bam <sample name> <output bam>", "Output bam file for the given sample"

      def ref(args)
        config[:reference_fa] = args.first
      end
      usage "ref <filename>", "Reference Fasta file to use for reordering."
    end
  end
end
