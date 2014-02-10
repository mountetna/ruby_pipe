module Genome
  class FixMate 
    include Pipeline::Step
    runs_tasks :sam_to_bam, :verify_mate, :enforce_label, :mark_duplicates
    runs_on :samples
    resources :walltime => 50

    class SamToBam 
      include Pipeline::Task
      requires_file :paired_merge_sam
      dumps_file :paired_merge_bam

      def run
        log_info "Converting SAM to BAM..."
        samtools "view -bS", config.paired_merge_sam, config.paired_merge_bam
      end
    end
 
    class VerifyMate 
      include Pipeline::Task
      requires_file :paired_merge_bam
      dumps_file :mated_bam

      def run
        log_info "Verifying mate information"
        picard :fix_mate_information, :INPUT => config.paired_merge_bam, :OUTPUT=> config.mated_bam or error_exit "Verify mate information failed"
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
                :CREATE_INDEX => :false, :SORT_ORDER => "coordinate",
                  :RGID => config.sample_name, :RGLB => config.sample_name,
                  :RGPL => config.platform, :RGPU => config.platform_unit, :RGSM => config.sample_name,
                  :VALIDATION_STRINGENCY => "LENIENT" or error_exit "Relabel failed"
      end
    end
    class MarkDuplicates
      include Pipeline::Task
      requires_file :raw_sample_bam
      outs_file :dedup_sample_bam, :duplication_metrics

      def run
      log_info "Mark duplicates"
      picard :mark_duplicates, :INPUT => config.raw_sample_bam,
              :OUTPUT => config.dedup_sample_bam, :METRICS_FILE => config.duplication_metrics, 
              :ASSUME_SORTED => :true,
              :REMOVE_DUPLICATES => :false, :CREATE_INDEX => :false or error_exit "Mark duplicates failed"
      end
    end
  end
end
