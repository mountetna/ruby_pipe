module Genome 
  class Merge
    include Pipeline::Step
    runs_tasks :merge_sam
    runs_on :samples#,:inputs
    resources :threads => 12, :walltime => 50

    class MergeSam
      include Pipeline::Task
      requires_files :aligned_bams
      dumps_file :raw_sample_bam

      def run
        log_info "Merging samll reads into one SAM file"
        picard :merge_sam_files, :O => config.raw_sample_bam, :I => config.aligned_bams, :CREATE_INDEX => :true or error_exit "BAM file merge failed."
      end
    end
  end

end
