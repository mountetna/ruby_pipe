module Genome 
  class SamMerge
    include Pipeline::Step
    runs_tasks :merge_sam
    runs_on :samples#,:inputs
    resources :threads => 12, :walltime => 50

    class MergeSam
      include Pipeline::Task
      requires_files :paired_sams
      dumps_file :paired_merge_sam

      def run
        log_info "Merging samll reads into one SAM file"
        picard :merge_sam_files, :CREATE_INDEX => :false, :SORT_ORDER => "unsorted", :USE_THREADING => "true",
          :OUTPUT => config.paired_merge_sam, :INPUT => config.paired_sams or error_exit "bam file merge failed."
      end
    end
  end

end
