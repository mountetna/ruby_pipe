module Genome 
  class Align 
    include Pipeline::Step
    runs_tasks :align_first, :align_second, :pair_reads
    runs_on :samples
    resources :threads => 12

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
  end
end
