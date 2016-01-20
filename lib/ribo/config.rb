module Ribo
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end

    def_var :output_bams do samples.map{|s| output_bam(s) } end

    dir_tree({
      ":tmp_dir" => {
        "@fraction_name.rsem.tmp" => :rsem_tmp_dir
      },
      ":scratch_dir" => {
        "@sample_name" => {
          "@fraction_name.clipped.fq" => :clipped_fastq,
          "@fraction_name.ribo.sam" => :ribo_sam,
          "@fraction_name.non_ribo.sam" => :non_ribo_sam,
          "@fraction_name.non_ribo.fq" => :non_ribo_fastq,
          "@fraction_name.read.sai" => :read_sai,
          "@fraction_name.mapped.sam" => :mapped_sam,
          "@fraction_name.genome.sorted.bam" => :genome_bam,
          "@fraction_name.unaligned.sam" => :unaligned_sam,
          "@fraction_name.unaligned.fq" => :unaligned_fastq,
          "@fraction_name.aligned.sam" => :bwa_aligned_sam,
          "@fraction_name.combined.header" => :combined_header,
          "@fraction_name.tophat.header" => :tophat_header,
          "@fraction_name.genome.header" => :genome_header,
          "@fraction_name.tophat.header" => :tophat_header,
          "@fraction_name.tophat" => {
            "." => :tophat_scratch,
            "accepted_hits.bam" => :accepted_bam,
            "unmapped.bam" => :unmapped_bam,
            "merged.bam" => :merged_bam,
            "tophat.bam" => :tophat_bam,
          },
          "@fraction_name.rsem" => {
            "." => :rsem_scratch_dir,
            "@fraction_name.genes.results" => :rsem_scratch_genes_results,
            "@fraction_name.isoforms.results" => :rsem_scratch_isoforms_results,
            "@fraction_name.transcript.bam" => :rsem_scratch_txp_unsorted_bam,
            "@fraction_name.transcript.sorted.bam" => :rsem_scratch_txp_bam,
            "@fraction_name.transcript.sorted.bam.bai" => :rsem_scratch_txp_bai,
            "@fraction_name.genome.sorted.bam" => :rsem_genome_bam,
            "@fraction_name.genome.sorted.bam.bai" => :rsem_genome_bai
          },
          "@fraction_name.rsem_aligned.bam" => :rsem_aligned_bam,
          "@fraction_name.rsem_aligned.sam" => :rsem_aligned_sam,
          "@fraction_name.rsem_fixed_aligned.bam" => :rsem_fixed_aligned_bam,
          "@fraction_name.rsem.header.sam" => :rsem_header,
          "@fraction_name.rsem.sorted.header.sam" => :sorted_rsem_header,
          "@fraction_name.bwa.header.sam" => :bwa_header,
          "@fraction_name.tophat_edit.bam" => :tophat_edit_bam,
          "@fraction_name.tophat_sort.bam" => :tophat_sort_bam,
          "@fraction_name.coverage.sam" => :coverage_sam
        }
      },
      ":output_dir" => {
        "transcript_model.gtf" => :transcript_model_gtf,
        "@sample_name" => {
          "@fraction_name.transcript_model.cov" => :transcript_model_coverage_base,
          "@fraction_name.normal.cov" => :normal_cov,
          "@fraction_name.null.cov" => :null_cov,
          "@fraction_name.bam" => :output_bam,
          "@fraction_name.bai" => :output_bai,
          "@fraction_name.aligned.bam" => :bwa_aligned_bam
        },
        "babel" => {
          "@babel_name" => :babel_output,
          "@babel_name.within.babel" => :within_babel,
          "@babel_name.combined.babel" => :combined_babel,
          "@babel_name.between.babel" => :between_babel
        },
        "@cohort_name.normal_cov" => :normal_summary,
        "@cohort_name.transcript_model.cov" => :transcript_model_summary_base,
        "@cohort_name.orf.cov" => :orf_summary,
        "@cohort_name.null_cov" => :null_summary,
      },
      ":metrics_dir" => {
        "@sample_name" => {
          "@fraction_name.flagstat" => :qc_flag,
          "@fraction_name.rnaseq_metrics" => :qc_rnaseq,
          "@fraction_name.rnaseq_pdf" => :qc_pdf,
          "@fraction_name.splice_counts" => :qc_splice_counts,
          "@fraction_name.alignment_metrics" => :qc_align_metrics,
          "@fraction_name.rrna_metrics" => :qc_rrna_metrics,
        },
        "@cohort_name.qc_summary" => :qc_summary
      }
    })

    job_items :babel

    def_var :fractions do ( samples.map(&:rna) + samples.map(&:rp) ).compact end

    # align
    def_var :input_fastq do job_item.input_fastq end

    def_var :babel_num_reps do 10000 end
    def_var :babel_min_rna do 10 end
    def_var :babel_min_rpkm do 0.2 end

    #coverage
    def_var :model_type do :unified_model end
    def_var :transcript_model_coverage do |name,fraction| 
      transcript_model_coverage_base(fraction)
        .sub(/transcript_model/,name.to_s)
    end
    def_var :transcript_model_coverages do 
      transcript_model_regions.map do |name|
        transcript_model_coverage(name)
      end
    end
    def_var :transcript_model_summary do |name| 
      transcript_model_summary_base.sub(/transcript_model/,name.to_s)
    end
    def_var :transcript_model_summaries do 
      transcript_model_regions.map do |name|
        transcript_model_summary(name)
      end
    end
    def_var :transcript_model_regions do [ :utr5, :utr3, :orf, :start, :stop ] end
    #def_var :transcript_model_regions do [ :orf ] end

    def init_hook
      samples.each do |s|
        s.rna.add_member :fraction_name, "#{s.sample_name}.rna" if s.rna
        s.rp.add_member :fraction_name, "#{s.sample_name}.rp" if s.rp
      end
    end
  end
end
