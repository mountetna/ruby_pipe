module Ribo
  class RsemAlign
    include Pipeline::Step
    runs_tasks :rsem_count, :collect_unmapped, :make_unaligned_fastq, :cull_aligned
    resources :threads => 12
    runs_on :fractions

    class RsemCount
      include Pipeline::Task
      requires_file :clipped_fastq
      outs_file :rsem_genome_bam
      
      def run
        rsem :calculate_expression,
            :temporary_folder => File.realdirpath(config.rsem_tmp_dir),
            :num_threads => config.threads,
            :output_genome_bam => true,
            :args => {
              :fq1s => rsem_format(config.clipped_fastq),
              :reference => rsem_format(config.reference_rsem),
              :sample_name => config.sample_replicate_name
            }, :output => config.rsem_scratch_dir or error_exit "Could not run RSEM"
      end
    end
    class CollectUnmapped
      include Pipeline::Task
      requires_file :rsem_genome_bam
      outs_file :unaligned_sam

      def run
        log_info "Culling unaligned reads"
        picard :view_sam,
          :I => config.rsem_genome_bam,
          :ALIGNMENT_STATUS => :Unaligned,
          :out => config.unaligned_sam or error_exit "picard view_sam failed"
      end
    end

    class MakeUnalignedFastq
      include Pipeline::Task
      requires_file :unaligned_sam
      dumps_file :unaligned_fastq

      def run
        picard :sam_to_fastq,
          :I => config.unaligned_sam,
          :FASTQ => config.unaligned_fastq or error_exit "picard sam_to_fastq failed"
      end
    end

    class CullAligned
      include Pipeline::Task
      requires_file :rsem_genome_bam
      outs_file  :rsem_aligned_bam

      def run
        log_info "Culling aligned reads"
        samtools "view -bh -F 4", config.rsem_genome_bam, config.rsem_aligned_bam
      end
    end
  end
  class BwaAlign
    include Pipeline::Step
    runs_on :fractions
    runs_tasks :align_unmapped, :sam_to_bam

    class AlignUnmapped
      include Pipeline::Task

      requires_file :unaligned_fastq
      dumps_file :bwa_aligned_sam

      def run
        bwa_mem fq1: config.unaligned_fastq, out: config.bwa_aligned_sam or error_exit "Could not run BWA"
      end
    end
    class SamToBam
      include Pipeline::Task

      requires_file :bwa_aligned_sam
      dumps_file :bwa_aligned_bam
      def run
        samtools "view -bh", config.bwa_aligned_sam, config.bwa_aligned_bam
      end
    end
  end
  class CombineRsem
    include Pipeline::Step
    runs_on :fractions
    runs_task :merge_reads

    class MergeReads
      include Pipeline::Task

      requires_files :bwa_aligned_bam, :rsem_aligned_bam
      outs_file :output_bam

      def run
        sam_merge config.output_bam, config.bwa_aligned_bam, config.rsem_aligned_bam
      end
    end
  end
end
