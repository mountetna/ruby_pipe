module Ribo
  class Coverage
    include Pipeline::Step
    runs_tasks :normal_coverage, :random_coverage
    job_list do config.samples end

    class NormalCoverage
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :normal_cov

      def run
        log_info "Mapping coverage to reference genes"
        samtools "view -h -q 1 -F 4", config.sample_bam, config.coverage_sam
        htseq_count :input => config.coverage_sam, :gtf => config.hg19_unified_gtf, :type => "unified_model", :out => config.normal_cov or error_exit "Computing normal coverage failed."
        File.unlink config.coverage_sam
        #coverage_bed config.sample_bam, config.reference_gtf, config.normal_cov or error_exit "Computing normal coverage failed."
      end
    end
    class RandomCoverage
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :random_cov

      def run
        log_info "Mapping coverage to randomized genes"
        coverage_bed config.sample_bam, config.hg19_null_gtf, config.random_cov or error_exit "Computing random coverage failed."
      end
    end
  end
end
