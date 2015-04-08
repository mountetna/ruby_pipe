module Ribo
  class Babel
    include Pipeline::Step
    runs_tasks :run_babel
    runs_on :babel_tests

    class NormalCoverage
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :normal_cov

      def run
        log_info "Mapping coverage to reference genes"
        samtools "view -h -q 1 -F 4", config.sample_bam, config.coverage_sam
        htseq_count :input => config.coverage_sam, :gtf => config.reference_unified_gtf, :type => config.model_type, :out => config.normal_cov or error_exit "Computing normal coverage failed."
        File.unlink config.coverage_sam
        #coverage_bed config.sample_bam, config.reference_gtf, config.normal_cov or error_exit "Computing normal coverage failed."
      end
    end
    class NullCoverage
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :null_cov

      def run
        log_info "Mapping coverage to randomized genes"
        coverage_bed config.sample_bam, config.reference_null_gtf, config.null_cov or error_exit "Computing random coverage failed."
      end
    end
  end
end
