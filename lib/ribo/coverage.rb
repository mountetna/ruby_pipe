module Ribo
  class Coverage
    include Pipeline::Step
    runs_tasks :normal_coverage, :transcript_model_coverage #, :null_coverage
    runs_on :fractions

    class NormalCoverage
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :normal_cov

      def run
        log_info "Mapping coverage to reference genes"
        run_cmd "samtools view -h -F 260 #{config.sample_bam} | awk '$5 != 0 || $0 ~ /NM:i:0/'  > #{config.coverage_sam}"
        htseq_count :input => config.coverage_sam, 
          :gtf => config.reference_unified_gtf, 
          :type => config.model_type, 
          :out => config.normal_cov or error_exit "Computing normal coverage failed."
        File.unlink config.coverage_sam
        #coverage_bed config.sample_bam, config.reference_gtf, config.normal_cov or error_exit "Computing normal coverage failed."
      end
    end
    class TranscriptModelCoverage
      include Pipeline::Task
      requires_file :transcript_model_gtf
      outs_files :transcript_model_coverages

      def run
        log_info "Mapping coverage to reference genes"
        run_cmd "samtools view -h -F 260 #{config.sample_bam} | awk '$5 != 0 || $0 ~ /NM:i:0/'  > #{config.coverage_sam}"

        config.transcript_model_regions.each do |region|
          htseq_count input: config.coverage_sam, 
            gtf: config.transcript_model_gtf, 
            type: region,
            out: config.transcript_model_coverage(region) or error_exit "Computing coverage for transcript region #{region} failed."
        end
        File.unlink config.coverage_sam
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
