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

    # align
    def_var :input_fastq do job_item.input_fastq end
    def_var :clipped_fastq do sample_scratch_file "clipped.fq" end
    def_var :read_sai do sample_scratch_file "read.sai" end
    def_var :mapped_sam do sample_scratch_file "mapped.sam" end
    def_var :genome_bam do sample_scratch_file "genome.bam" end

    # tophat
    def_var :unaligned_sam do sample_scratch_file "unaligned.sam" end
    def_var :unaligned_fastq do sample_scratch_file "unaligned.fq" end
    def_var :aligned_bam  do sample_scratch_file "aligned.bam" end
    def_var :combined_header do sample_scratch_file "combined.header" end
    def_var :tophat_header do sample_scratch_file "tophat.header" end
    def_var :genome_header do sample_scratch_file "genome.header" end
    def_var :tophat_scratch do sample_scratch_file "tophat" end
    def_var :tophat_file do |aff| File.join tophat_scratch, aff end
    def_var :accepted_bam do tophat_file "accepted_hits.bam" end
    def_var :unmapped_bam do tophat_file "unmapped.bam" end
    def_var :merged_bam do tophat_file "merged.bam" end
    def_var :tophat_bam do sample_scratch_file "tophat.bam" end

    # qc
    def_var :qc_flag do |s| sample_metrics_file "flagstat" end
    def_var :qc_rnaseq do |s| sample_metrics_file "rnaseq_metrics" end
    def_var :qc_pdf do |s| sample_metrics_file "rnaseq_pdf" end
    def_var :qc_align_metrics do |s| sample_metrics_file "alignment_metrics" end

    # combine
    def_var :tophat_edit_bam do sample_scratch_file "tophat_edit.bam" end
    def_var :tophat_sort_bam do sample_scratch_file "tophat_sort.bam" end

    #coverage
    def_var :normal_cov do sample_output_file "normal.cov" end
    def_var :random_cov do sample_output_file "null.cov" end
    def_var :coverage_sam do sample_scratch_file "coverage.sam" end
  end
end
