module Genome
  class FixMate 
    include Pipeline::Step
    runs_tasks :verify_mate, :enforce_label
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
  end
end
