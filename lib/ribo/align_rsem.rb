require 'sam'

module Ribo
  class RsemAlign
    include Pipeline::Step
    runs_tasks :rsem_count, :collect_unmapped, :make_unaligned_fastq, :copy_rsem_bam
    resources :threads => 12, :memory => "4gb"
    runs_on :fractions

    class RsemCount
      include Pipeline::Task
      requires_file :non_ribo_fastq
      outs_file :rsem_genome_bam
      
      def run
        rsem :calculate_expression,
            :temporary_folder => config.rsem_tmp_dir,
            :num_threads => config.threads,
            :output_genome_bam => true,
            :strand_specific => true,
            :bowtie2 => true,
            :args => {
              :fq1s => rsem_format(config.non_ribo_fastq),
              :reference => rsem_format(config.reference_rsem),
              :sample_name => config.fraction_name,
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

    class CopyRsemBam
      include Pipeline::Task
      requires_file :rsem_genome_bam, :rsem_genome_bai
      outs_file :output_bam, :output_bai

      def run
        begin
          FileUtils.cp config.rsem_genome_bam, config.output_bam
        rescue
          error_exit "Failed to copy BAM file"
        end

        begin
          FileUtils.cp config.rsem_genome_bai, config.output_bai
        rescue
          error_exit "Failed to copy BAI file"
        end
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
      outs_file :bwa_aligned_bam
      def run
        samtools "view -bSh", config.bwa_aligned_sam, config.bwa_aligned_bam
      end
    end
  end
end
