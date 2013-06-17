require 'pipeline'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    dir_tree({
      ":output_dir" => {
        "@sample_name.reordered.bai" => :sample_bai
      }
    })
  end

  class IndexBam
    include Pipeline::Step
    runs_tasks :index_bam
    runs_on :samples

    class IndexBam
      include Pipeline::Task
      requires_file :input_bam
      outs_file :sample_bai

      def run
        log_info "Indexing bam file"
        samtools "index", config.input_bam or error_exit "Couldn't index bam"
      end
    end
  end
  class Index
    include Pipeline::Script
    runs_steps :index_bam
    def_module :default, :index_bam => true
    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input_bam(args)
        sample = get_sample args.first
        sample[:input_bam] = args[1]
      end
      usage "input_bam <sample name> <input_bam>", "Input bam file for the given sample"
    end
  end
end
