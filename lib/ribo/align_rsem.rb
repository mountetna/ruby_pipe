require 'sam'

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
            :bowtie2 => true,
            :args => {
              :fq1s => rsem_format(config.clipped_fastq),
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
        samtools "view -bSh", config.bwa_aligned_sam, config.bwa_aligned_bam
      end
    end
  end
  class CombineRsem
    include Pipeline::Step
    runs_on :fractions
    runs_task :build_header, :build_reads, :build_bam, :merge_aligned_bams

    class BuildHeader
      include Pipeline::Task
      requires_file :bwa_aligned_bam
      dumps_file :sorted_rsem_header

      def run
        samtools "view -H", config.bwa_aligned_bam, config.bwa_header
        samtools "view -H", config.rsem_aligned_bam, config.rsem_header

        bwa_sam = Sam.new.parse_sam config.bwa_header
        rsem_sam = Sam.new.parse_sam config.rsem_header

        rsem_sam.header.sequences.sort_by! do |rseq|
          bwa_sam.header.sequences.index do |bseq|
            bseq.sn == rseq.sn
          end
        end

        rsem_sam.print config.sorted_rsem_header
        File.unlink config.rsem_header, config.bwa_header
      end
    end

    class BuildReads
      include Pipeline::Task

      requires_file :rsem_aligned_bam
      dumps_file :rsem_aligned_sam

      def run
        samtools :view, config.rsem_aligned_bam, config.rsem_aligned_sam
      end
    end

    class BuildBam
      include Pipeline::Task

      requires_files :sorted_rsem_header, :rsem_aligned_sam
      dumps_file :rsem_fixed_aligned_bam

      def run
        run_cmd "cat #{config.sorted_rsem_header} #{config.rsem_aligned_sam} | samtools view -bSh - > #{config.rsem_fixed_aligned_bam}" or error_exit "Could not merge reads and header"
      end
    end

    class MergeAlignedBams
      include Pipeline::Task

      requires_files :bwa_aligned_bam, :rsem_fixed_aligned_bam
      outs_file :output_bam

      def run
        log_info "Merging unmapped and accepted reads"
        picard :merge_sam_files, :CREATE_INDEX => :true, :MERGE_SEQUENCE_DICTIONARIES => :true, :SORT_ORDER => :coordinate, :O => config.output_bam, :I => [ config.bwa_aligned_bam, config.rsem_fixed_aligned_bam ] or error_exit "Could not merge bam."
      end
    end
  end
end
