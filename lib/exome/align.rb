module Exome
  class Align 
    include Pipeline::Step
    runs_tasks :align_first, :align_second, :pair_reads, :verify_mate, :enforce_label
    job_list do config.samples.collect(&:inputs).flatten end
    resources :threads => 6

    class AlignFirst 
      include Pipeline::Task
      requires_file :input_fastq1
      dumps_file :read1_sai

      def run
        log_info "Aligning first-in-pair reads"
        bwa_aln :fq => config.input_fastq1, :out => config.read1_sai or error_exit "First BWA alignment failed"
      end
    end

    class AlignSecond 
      include Pipeline::Task
      requires_file :input_fastq2
      dumps_file :read2_sai

      def run
        log_info "Aligning second-in-pair reads"
        bwa_aln :fq => config.input_fastq2, :out => config.read2_sai or error_exit "Second BWA alignment failed"
      end
    end

    class PairReads 
      include Pipeline::Task
      requires_files :read1_sai, :read2_sai
      dumps_file :paired_sam

      def run
        log_info "Pairing aligned reads"
        bwa_pair :m1 => config.read1_sai, :m2 => config.read2_sai, 
          :fq1 =>  config.input_fastq1, :fq2 => config.input_fastq2, :out => config.paired_sam or error_exit "BWA sampe failed"
      end
    end

    class VerifyMate 
      include Pipeline::Task
      requires_file :paired_sam
      dumps_file :mated_sam

      def run
        log_info "Verifying mate information"
        picard :fix_mate_information, :INPUT => config.paired_sam, :OUTPUT=> config.mated_sam or error_exit "Verify mate information failed"
      end
    end

    class EnforceLabel 
      include Pipeline::Task
      requires_file :mated_sam
      dumps_file :aligned_bam

      def run
        log_info "Enforce read group assignments"
        picard :add_or_replace_read_groups, :INPUT => config.mated_sam,
                :OUTPUT => config.aligned_bam,
                :CREATE_INDEX => :true,
                  :RGID => config.sample_name, :RGLB => config.sample_name,
                  :RGPL => config.platform, :RGPU => config.platform_unit, :RGSM => config.sample_name,
                  :VALIDATION_STRINGENCY => "LENIENT" or error_exit "Relabel failed"
      end
    end
  end
end
