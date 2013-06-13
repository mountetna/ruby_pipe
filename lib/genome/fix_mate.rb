module Genome
  class FixMate 
    include Pipeline::Step
    runs_tasks :verify_mate, :enforce_label, :mark_duplicates 
    runs_on :samples

    class VerifyMate 
      include Pipeline::Task
      requires_file :input_bam
      dumps_file :mated_bam

      def run
        log_info "Verifying mate information"
        picard :fix_mate_information, :INPUT => config.input_bam, :OUTPUT=> config.mated_bam or error_exit "Verify mate information failed"
      end
    end

    class EnforceLabel 
      include Pipeline::Task
      requires_file :mated_bam
      outs_file :raw_sample_bam

      def run
        log_info "Enforce read group assignments"
        picard :add_or_replace_read_groups, :INPUT => config.mated_bam,
                :OUTPUT => config.raw_sample_bam,
                :CREATE_INDEX => :true,
                  :RGID => config.sample_name, :RGLB => config.sample_name,
                  :RGPL => config.platform, :RGPU => config.platform_unit, :RGSM => config.sample_name,
                  :VALIDATION_STRINGENCY => "LENIENT" or error_exit "Relabel failed"
      end
    end
    class MarkDuplicates
      include Pipeline::Task
      requires_file :raw_sample_bam #:merged_library_bam
      outs_file :dedup_sample_bam, :duplication_metrics

      def run
      log_info "Mark duplicates"
      picard :mark_duplicates, :INPUT => config.raw_sample_bam,#config.merged_library_bam,
              :OUTPUT => config.dedup_sample_bam, :METRICS_FILE => config.duplication_metrics, 
              :REMOVE_DUPLICATES => :false, :CREATE_INDEX => :true or error_exit "Mark duplicates failed"
      end
    end
  end
end
