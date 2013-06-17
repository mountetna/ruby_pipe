require 'pipeline'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def_var :feature_name do cohort_name end

    dir_tree({
      ":output_dir" => {
        "@sample_name.:feature_name.cov" => :feature_cov
      }
    })
  end

  class ComputeCoverage
    include Pipeline::Step
    runs_tasks :gtf_coverage
    runs_on :samples

    class GtfCoverage
      include Pipeline::Task
      requires_file :input_bam
      outs_file :feature_cov

      def run
        log_info "Mapping coverage to randomized genes"
        coverage_bed config.input_bam, config.reference_gtf, config.feature_cov or error_exit "Computing random coverage failed."
      end
    end
  end
  class Coverage
    include Pipeline::Script
    runs_steps :compute_coverage
    def_module :default, :compute_coverage => true
    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.first
        sample[:input_bam] = args[1]
      end
      usage "input <sample name> <input_bam>", "Input bam file for the given sample"

      def feature_name(args)
        config[:feature_name] = args.first
      end
      usage "feature_name <name>", "Specify the feature name for output files."

      def gtf(args)
        config[:reference_gtf] = args.first
      end
      usage "gtf <filename>", "Reference GTF file to use for coverage features."
    end
  end
end
