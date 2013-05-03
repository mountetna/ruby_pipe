module Ribo
  class Combine
    include Pipeline::Step
    runs_tasks :sort_seq, :enforce_headers, :merge_bams
    resources :threads => 12
    runs_on :samples

    class SortSeq
      include Pipeline::Task
      requires_files :tophat_bam
      dumps_files :tophat_sort_bam

      def run
        log_info "Sorting sequence dictionary"
        picard :reorder_sam, :I => config.tophat_bam, :O => config.tophat_sort_bam, :REFERENCE => config.hg19_fa or error_exit "Could not reorder bam."
      end
    end

    class EnforceHeaders
      include Pipeline::Task
      requires_file :tophat_sort_bam
      dumps_file :tophat_edit_bam

      def run
        samtools "view -H", config.genome_bam, config.genome_header or error_exit "samtools view header failed"
        samtools "view -H", config.tophat_sort_bam, config.tophat_header or error_exit "samtools view header failed"
        File.open(config.combined_header,"w") do |f|
          f.puts File.read(config.genome_header)
          f.puts File.foreach(config.tophat_header).grep(/ID:TopHat/)
        end

        picard :replace_sam_header, :HEADER => config.combined_header,
          :INPUT => config.tophat_sort_bam,
          :OUTPUT => config.tophat_edit_bam or error_exit "picard replace_sam_header failed"

        File.unlink config.genome_header
        File.unlink config.tophat_header
        File.unlink config.combined_header
      end
    end


    class MergeBams
      include Pipeline::Task
      requires_file :tophat_edit_bam, :aligned_bam
      outs_file :output_bam

      def run
        log_info "Merge genome and spliced alignments"
        picard :merge_sam_files,
          :SO => :coordinate,
          :OUTPUT => config.output_bam,
          :INPUT => [ config.aligned_bam, config.tophat_edit_bam ],
          :MSD => :true or error_exit "picard merge failed"

        sam_index config.output_bam or error_exit "index failed"
      end
    end
  end
end

