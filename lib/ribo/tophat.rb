module Ribo
  class Tophat
    include Pipeline::Step
    runs_tasks :cull_unaligned, :cull_aligned, :make_unaligned_fastq, :tophat_align, :merge_reads, :enforce_label
    resources :threads => 12
    runs_on :samples

    class CullUnaligned
      include Pipeline::Task
      requires_file :genome_bam
      dumps_file :unaligned_sam

      def run
        log_info "Culling unaligned reads"
        picard :view_sam,
          :I => config.genome_bam,
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
      requires_file :genome_bam
      outs_file  :aligned_bam

      def run
        log_info "Culling aligned reads"
        samtools "view -bh -F 4", config.genome_bam, config.aligned_bam
      end
    end

    class TophatAlign
      include Pipeline::Task
      requires_files :unaligned_fastq
      dumps_files :accepted_bam, :unmapped_bam

      def run
        log_info "Running tophat"
        tophat :max_multihits => 1,
          :output_dir => config.tophat_scratch,
          :min_anchor =>  6,
          :splice_mismatches => 1,
          :no_novel_juncs => true,
          :library_type => "fr-unstranded",
          :transcriptome_index => config.transcriptome_idx,
          :fq1 => config.unaligned_fastq or error_exit "tophat failed"
      end
    end

    class MergeReads
      include Pipeline::Task
      requires_files :accepted_bam, :unmapped_bam
      dumps_files :merged_bam

      def run
        log_info "Merging unmapped and accepted reads"
        picard :merge_sam_files, :O => config.merged_bam, :I => [ config.accepted_bam, config.unmapped_bam ], :MSD => :true or error_exit "Could not merge bam."
      end
    end

    class EnforceLabel 
      include Pipeline::Task
      requires_file :merged_bam
      outs_file :tophat_bam

      def run
        log_info "Enforce read group assignments"
        picard :add_or_replace_read_groups, :INPUT => config.merged_bam,
                :CREATE_INDEX => :true,
                :SO => :coordinate,
                :OUTPUT => config.tophat_bam,
                  :RGID => config.sample_name, :RGLB => config.sample_name,
                  :RGPL => :illumina, :RGPU => :None, :RGSM => config.sample_name,
                  :CN => "Babel@UCSF" or error_exit "Relabel failed"
      end
    end
  end
end

