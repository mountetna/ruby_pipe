module Genome
  class DumpFastqs
    include Pipeline::Step
    runs_task :purify_bam, :create_fastq_from_bam
    runs_on :samples, :inputs

    class PurifyBam
      include Pipeline::Task
      skip_without :reads_bam
      outs_file :chaste_bam

      def run
        samtools "view -bh -F 0x200", config.reads_bam, config.chaste_bam or error_exit "Couldn't filter BAM"
      end
    end

    class CreateFastqFromBam
      # you only need to do this if you have your bams in reads
      include Pipeline::Task
      skip_without :reads_bam
      requires_file :chaste_bam
      outs_file :reads1_fastq, :reads2_fastq

      def run
        picard :sam_to_fastq, :INPUT => config.chaste_bam, :FASTQ => ">(gzip -c > #{config.reads1_fastq})", :SECOND_END_FASTQ => ">(gzip -c > #{config.reads2_fastq})" or error_exit "Could not split bam file to fastqs"
      end
    end
  end

  class CombineFastqs
    include Pipeline::Step
    runs_task :compute_chunks
    runs_on :samples, :inputs

    class ComputeChunks
      include Pipeline::Task
      requires_files :input_fastq1, :input_fastq2
      outs_file :chunk_info

      def run
        lines1 =  read_cmd("wc -l <(zcat #{config.input_fastq1})") or error_exit "Could not count fastq files"
        lines1 = lines1.split.first.to_i 
        lines2 =  read_cmd("wc -l <(zcat #{config.input_fastq2})") or error_exit "Could not count fastq files"
        lines2 = lines2.split.first.to_i 
        error_exit "Fastq files differ in size" if lines1 != lines2
        error_exit "Fastq files are not properly formatted" if lines1 % 4 != 0
        File.open(config.chunk_info,"w"){|f| f.puts lines1/4}
      end
    end
  end

  class Align 
    include Pipeline::Step
    runs_tasks :make_fastq_chunk, :align_first, :align_second, :pair_reads, :verify_mate, :mark_duplicates, :enforce_label
    has_tasks :make_fastq_chunk, :align_first, :align_second, :pair_reads, :align_mem, :verify_mate, :mark_duplicates, :enforce_label
    runs_on :samples, :inputs, :chunks
    resources :threads => 12

    class MakeFastqChunk
      include Pipeline::Task
      requires_files :input_fastq1, :input_fastq2
      dumps_file :chunk1_fastq, :chunk2_fastq

      def run
        # get some lines from the middle
        chunk_start = config.chunk_size * config.chunk.chunk_number + 1
        chunk_end = chunk_start + config.chunk_size - 1
        run_cmd "sed -n '#{chunk_start}, #{chunk_end}p; #{chunk_end}q' <(zcat #{config.input_fastq1}) | gzip -c > #{config.chunk1_fastq}" or error_exit "Couldn't chunk file"
        run_cmd "sed -n '#{chunk_start}, #{chunk_end}p; #{chunk_end}q' <(zcat #{config.input_fastq2}) | gzip -c > #{config.chunk2_fastq}" or error_exit "Couldn't chunk file"
      end
    end

    class AlignFirst 
      include Pipeline::Task
      requires_file :chunk1_fastq
      dumps_file :read1_sai

      def run
        log_info "Aligning first-in-pair reads"
        bwa_aln :fq => config.chunk1_fastq, :out => config.read1_sai or error_exit "First BWA alignment failed"
      end
    end

    class AlignSecond 
      include Pipeline::Task
      requires_file :chunk2_fastq
      dumps_file :read2_sai

      def run
        log_info "Aligning second-in-pair reads"
        bwa_aln :fq => config.chunk2_fastq, :out => config.read2_sai or error_exit "Second BWA alignment failed"
      end
    end

    class PairReads 
      include Pipeline::Task
      requires_files :read1_sai, :read2_sai
      dumps_file :paired_sam

      def run
        log_info "Pairing aligned reads"
        bwa_pair :m1 => config.read1_sai, :m2 => config.read2_sai, 
          :fq1 =>  config.chunk1_fastq, :fq2 => config.chunk2_fastq, :out => config.paired_sam or error_exit "BWA sampe failed"
      end
    end

    class AlignMem 
      include Pipeline::Task
      requires_files :chunk1_fastq, :chunk2_fastq
      dumps_file :paired_sam

      def run
        log_info "Pairing aligned reads"
        bwa_mem :fq1 =>  config.chunk1_fastq, :fq2 => config.chunk2_fastq, :out => config.paired_sam or error_exit "BWA mem failed"
      end
    end

    class VerifyMate 
      include Pipeline::Task
      requires_file :paired_sam
      dumps_file :mated_bam

      def run
        log_info "Verifying mate information"
        picard :fix_mate_information, :INPUT => config.paired_sam, 
                :SORT_ORDER => :coordinate,
                :OUTPUT=> config.mated_bam or error_exit "Verify mate information failed"
      end
    end

    class MarkDuplicates
      include Pipeline::Task
      requires_file :mated_bam
      outs_file :dedup_bam, :duplication_metrics

      def run
	log_info "Mark duplicates"
	picard :mark_duplicates, :INPUT => config.mated_bam,
          :OUTPUT => config.dedup_bam, :METRICS_FILE => config.duplication_metrics, 
          :REMOVE_DUPLICATES => :false or error_exit "Mark duplicates failed"
      end
    end

    class EnforceLabel 
      include Pipeline::Task
      requires_file :dedup_bam
      outs_file :aligned_bam

      def run
        log_info "Enforce read group assignments"
        picard :add_or_replace_read_groups, :INPUT => config.dedup_bam,
                :OUTPUT => config.aligned_bam,
                :SORT_ORDER => :coordinate,
                :CREATE_INDEX => :true,
                  :RGID => config.sample_name, :RGLB => config.sample_name,
                  :RGPL => config.platform, :RGPU => config.platform_unit, :RGSM => config.sample_name,
                  :VALIDATION_STRINGENCY => "LENIENT" or error_exit "Relabel failed"
      end
    end
  end
end
