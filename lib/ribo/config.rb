module Ribo
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end

    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :input_bam do |s| s ? sample(s)[:input_bam] : nil end
    def_var :output_bam do sample_output_file "#{sample_name}.bam" end
    def_var :output_bams do sample_names.map{|s| output_bam(s) } end

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
        "normal.cov" => :normal_cov,
        "null.cov" => :random_cov,
        "coverage.sam" => :coverage_sam,
      },
      ":metrics_dir" => {
        "@sample_name" => {
          "@sample_name.flagstat" => :qc_flag,
          "@sample_name.rnaseq_metrics" => :qc_rnaseq,
          "@sample_name.rnaseq_pdf" => :qc_pdf,
          "@sample_name.alignment_metrics" => :qc_align_metrics,
        }
      }
    })

    # align
    def_var :input_fastq do job_item.input_fastq end

    #coverage
  end
end
