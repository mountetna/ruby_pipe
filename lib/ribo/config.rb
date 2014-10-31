module Ribo
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end

    def_var :output_bams do samples.map{|s| output_bam(s) } end

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "clipped.fq" => :clipped_fastq,
          "read.sai" => :read_sai,
          "mapped.sam" => :mapped_sam,
          "genome.bam" => :genome_bam,
          "unaligned.sam" => :unaligned_sam,
          "unaligned.fq" => :unaligned_fastq,
          "aligned.bam" => :aligned_bam,
          "combined.header" => :combined_header,
          "tophat.header" => :tophat_header,
          "genome.header" => :genome_header,
          "tophat" => {
            "." => :tophat_scratch,
            "accepted_hits.bam" => :accepted_bam,
            "unmapped.bam" => :unmapped_bam,
            "merged.bam" => :merged_bam,
            "tophat.bam" => :tophat_bam,
          },
          "tophat_edit.bam" => :tophat_edit_bam,
          "tophat_sort.bam" => :tophat_sort_bam,
        }
      },
      ":output_dir" => {
        "@sample_name" => {
          "normal.cov" => :normal_cov,
          "null.cov" => :null_cov,
          "coverage.sam" => :coverage_sam,
          "@sample_name.bam" => :output_bam
        },
        "@cohort_name.normal_cov" => :normal_summary,
        "@cohort_name.null_cov" => :null_summary,
      },
      ":metrics_dir" => {
        "@sample_name" => {
          "@sample_name.flagstat" => :qc_flag,
          "@sample_name.rnaseq_metrics" => :qc_rnaseq,
          "@sample_name.rnaseq_pdf" => :qc_pdf,
          "@sample_name.splice_counts" => :qc_splice_counts,
          "@sample_name.alignment_metrics" => :qc_align_metrics,
        },
        "@cohort_name.qc_summary" => :qc_summary
      }
    })

    # align
    def_var :input_fastq do job_item.input_fastq end

    #coverage
    def_var :model_type do :unified_model end
  end
end
