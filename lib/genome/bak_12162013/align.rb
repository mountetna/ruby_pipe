module Genome 
  class CreateFastq
    include Pipeline::Step
    runs_tasks :create_fastq 
    runs_on :samples

    class CreateFastq
      include Pipeline::Task
      requires_file :pre_sample_bam
      outs_file :reads1_fastq, :reads2_fastq

      def run
        log_info "Converting BAM to FastQs (reads 1 and 2)"
        picard :sam_to_fastq, :INPUT => config.chaste_bam, :FASTQ => ">(gzip -c > #{config.reads1_fastq})", :SECOND_END_FASTQ => ">(gzip -c > #{config.reads2_fastq})" or error_exit "Could not split bam file to fastqs"
      end
    end
  end 

  class Align 
    include Pipeline::Step
    runs_tasks :align_first, :align_second, :pair_reads#, :verify_mate
    runs_on :samples, :inputs#:chunks
    resources :threads => 12
  
    class ChunkFastq
      include Pipeline::Task
      requires_files :reads1_fastq, :reads2_fastq
      dumps_file :chunk1_fastq, :chunk2_fastq
     
      def run
        log_info "Chunking fastqs...chunk chunk chunk"
        run_cmd "split -l 5000000 -a 3 -d <(zcat #{config.reads1_fastqs.join(" ")}) | gzip -c > #{config.chunk1_fastq}" or error_exit "Couldn't chunk file"
        run_cmd "split -l 5000000 -a 3 -d <(zcat #{config.reads2_fastqs.join(" ")}) | gzip -c > #{config.chunk2_fastq}" or error_exit "Couldn't chunk file"
      end
    end

    class AlignFirst 
      include Pipeline::Task
      requires_file :input_fastq1 #:chunk1_fastq
      dumps_file :read1_sai

      def run
        log_info "Aligning first-in-pair reads"
        bwa_aln :fq => config.input_fastq1, :out => config.read1_sai or error_exit "First BWA alignment failed"
      end
    end

    class AlignSecond 
      include Pipeline::Task
      requires_file :input_fastq2 #:chunk2_fastq
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

#    class VerifyMate 
#      include Pipeline::Task
#      requires_file :
#      dumps_file :mated_sam

#      def run
#        log_info "Verifying mate information"
#        picard :fix_mate_information, :INPUT => config.paired_sam, :OUTPUT=> config.mated_sam or error_exit "Verify mate information failed"
#      end
#    end

  end
end
