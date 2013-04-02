#!/usr/bin/env ruby
module Rna
  class TophatAlign
    include Pipeline::Step
    runs_tasks :tophat_align, :merge_reads, :sort_seq, :relabel_bam
    resources :threads => 12
    job_list do config.replicates end

    class TophatAlign
      include Pipeline::Task
      requires_files :input_fastq1s, :input_fastq2s
      dumps_files :accepted_bam, :unmapped_bam

      def run
        log_info "Running tophat"
        tophat :output_dir => config.tophat_scratch,
          :GTF => config.reference_gtf,
          :fq1 => config.input_fastq1s.join(","),
          :fq2 => config.input_fastq2s.join(","),
          :mate_inner_dist => config.frag_size or error_exit "Tophat failed."
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
    class SortSeq
      include Pipeline::Task
      requires_files :merged_bam
      dumps_files :sorted_bam

      def run
        log_info "Sorting sequence dictionary"
        #replace_dict config.merged_bam, config.hg19_dict, config.sorted_header or error_exit "Could not replace dictionary."
        #samtools "reheader #{config.sorted_header}", config.merged_bam
        picard :reorder_sam, :I => config.merged_bam, :O => config.sorted_bam, :REFERENCE => config.hg19_fa or error_exit "Could not reorder bam."
      end
    end
    class RelabelBam
      include Pipeline::Task
      requires_files :sorted_bam
      outs_file :output_bam

      def run
        log_info "Relabeling and sorting alignments"
        picard :add_or_replace_read_groups, :CREATE_INDEX => :true, :I => config.sorted_bam, :O => config.output_bam, :ID => config.sample_name, :LB => config.sample_name, :PL => :Illumina, :PU => :None,
          :SM => config.sample_name, :CN => :"Babel@UCSF" or error_exit "Could not label groups."
      end
    end
  end
end
