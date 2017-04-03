module Rna
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def init_hook
      # add various bells and whistles here
      samples.each do |s|
        s.replicates.each do |r|
          r.add_member :replicate_name, "r#{r.index}" if !r.replicate_name
        end
      end
    end

    empty_var :splice_gtf

    def_var :bam_label do "aligned.deduplicated.sorted" end
    def_var :fdr_cutoff do 0.10 end
    def_var :normal_bams do replicate_bams normal end

    def_var :replicate_bams do |s| (s||sample).replicates.map{ |r| replicate_bam r} end
    def_var :replicate_bam do |r| input_bam(r) || output_bam(r) end

    def_var :replicate_name do |r| (r || job_item).property :replicate_name end
    def_var :replicates do samples.collect(&:replicates).flatten end

    def_var :diff_exps do samples.collect(&:diff_exp).compact.flatten end

    def_var :sample_replicate_name do |r| "#{(r || job_item).property :sample_name}.#{(r || job_item).property :replicate_name}" end

    dir_tree({
      ":tmp_dir" => {
        "@cohort_name.@sample_name.@replicate_name.rsem_tmp" => :rsem_tmp_dir
        "@cohort_name.@sample_name.@replicate_name.rrna_tmp" => :rrna_tmp
      },
      ":scratch_dir" => {
        "@sample_name" => {
          "@replicate_name" => {
            "@sample_name.@replicate_name.rrna.sam" => :rrna_sam,
            "@sample_name.@replicate_name.rrna.bam" => :rrna_bam,
            "@sample_name.@replicate_name.non_rrna.sam" => :non_rrna_sam,
            "@sample_name.@replicate_name.non_rrna1.fq" => :non_rrna1_fastq,
            "@sample_name.@replicate_name.non_rrna2.fq" => :non_rrna2_fastq,
            "@sample_name.@replicate_name.non_rrna1.fastq.gz" => :non_rrna1_fastq_gz,
            "@sample_name.@replicate_name.non_rrna2.fastq.gz" => :non_rrna2_fastq_gz,
            "coverage.sam" => :coverage_sam,
            "splice.sam" => :splice_sam,
            "chimera1.fastq" => :chimera_fq1,
            "chimera2.fastq" => :chimera_fq2,
            "chimerascan" => {
              "." => :chimera_dir,
              "chimeras.bedpe" => :chimera_bedpe
            },
            "kallisto" => {
              "." => :kallisto_scratch,
              "abundance.tsv" => :kallisto_abundance_tsv,
            },
            "tophat" => {
              "." => :tophat_scratch,
              "accepted_hits.bam" => :accepted_bam,
              "unmapped.bam" => :unmapped_bam,
              "merged.bam" => :merged_bam,
              "merged.sorted.bam" => :sorted_bam,
              "sorted_header.txt" => :sorted_header,
            },
            "cufflinks" => {
              "." => :cufflinks_scratch,
              "transcripts.gtf" => :transcripts_gtf,
              "genes.fpkm_tracking" => :gene_tracking
            },
            "rsem" => {
              "." => :rsem_scratch_dir,
              #"tmp" => :rsem_tmp_dir,
              "@sample_name.@replicate_name.genes.results" => :rsem_scratch_genes_results,
              "@sample_name.@replicate_name.isoforms.results" => :rsem_scratch_isoforms_results,
              "@sample_name.@replicate_name.transcript.bam" => :rsem_scratch_txp_unsorted_bam,
              "@sample_name.@replicate_name.transcript.sorted.bam" => :rsem_scratch_txp_bam,
              "@sample_name.@replicate_name.transcript.sorted.bam.bai" => :rsem_scratch_txp_bai,
              "@sample_name.@replicate_name.genome.bam" => :rsem_scratch_genome_unsorted_bam,
              "@sample_name.@replicate_name.genome.header" => :rsem_scratch_genome_header,
              "@sample_name.@replicate_name.genome.sorted.bam" => :rsem_scratch_genome_sorted_bam,
              "@sample_name.@replicate_name.genome.patched.bam" => :rsem_scratch_genome_patched_bam,
              "@sample_name.@replicate_name.genome.sorted.bam.bai" => :rsem_scratch_genome_bai,
            },
            "@sample_name.@replicate_name.rsem.genome.sorted.bam" => :rsem_genome_sorted_bam,
            "@sample_name.@replicate_name.unaligned.bam" => :unaligned_bam,
            "@sample_name.@replicate_name.unaligned1.fq" => :unaligned1_fastq,
            "@sample_name.@replicate_name.unaligned2.fq" => :unaligned2_fastq,
            "@sample_name.@replicate_name.unaligned1.fastq.gz" => :unaligned1_fastq_gz,
            "@sample_name.@replicate_name.unaligned2.fastq.gz" => :unaligned2_fastq_gz,
            "@sample_name.@replicate_name.bwa_aligned.sam" => :bwa_aligned_sam
          },
          "cuffdiff_@normal_name" => {
            "." => :cuffdiff_dir,
            "gene_exp.diff" => :gene_exp_diff
          }
        },
        "@cohort_name" => {
          "cuffcompare" => {
            "cuffcomp" => :cuffcompare_prefix,
            "cuffcomp.tracking" => :tracking_file
          },
          "assembly.list" => :assembly_list,
          "merged.gtf" => :assembly_gtf
        }
      },
      ":metrics_dir" => {
        "@sample_name" => {
          "@sample_name.@replicate_name.flagstat" => :qc_flag,
          "@sample_name.@replicate_name.rnaseq_metrics" => :qc_rnaseq,
          "@sample_name.@replicate_name.rnaseq.pdf" => :qc_pdf,
          "@sample_name.@replicate_name.alignment_metrics" => :qc_align_metrics,
          "@sample_name.@replicate_name.bwa_flagstat" => :bwa_qc_flag,
          "@sample_name.@replicate_name.bwa_rnaseq_metrics" => :bwa_qc_rnaseq,
          "@sample_name.@replicate_name.bwa_rnaseq.pdf" => :bwa_qc_pdf,
          "@sample_name.@replicate_name.duplication_metrics" => :duplication_metrics,
          "@sample_name.@replicate_name.exon_total_cov" => :exon_total_cov,
          "@sample_name.@replicate_name.exon_unique_cov" => :exon_unique_cov,
          "@sample_name.@replicate_name.rrna_metrics" => :qc_rrna_metrics,
        },
        "@cohort_name.qc_summary" => :qc_summary
      },
      ":output_dir" => {
        "@sample_name" => {
          "@sample_name.@replicate_name.transcripts.gtf" => :output_gtf,
          "@sample_name.@replicate_name.:bam_label.bam" => :output_bam,
          "@sample_name.@replicate_name.:bam_label.unpatched.bam" => :output_unpatched_bam,
          "@sample_name.@replicate_name.header" => :output_header,
          "@sample_name.@replicate_name.:bam_label.bai" => :output_bai,
          "@sample_name.@replicate_name.genes.results" => :rsem_genes_results,
          "@sample_name.@replicate_name.isoforms.results" => :rsem_isoforms_results,
          "@sample_name.mutations" => :sample_mutations,
          "@sample_name.@normal_name.diff_exp" => :diff_exp_table,
          "@sample_name.@replicate_name.transcripts.cov" => :transcripts_cov,
          "@sample_name.@replicate_name.exon_splice_counts" => :exon_splice_counts,
          "@sample_name.@replicate_name.kallisto.genes.results" => :kallisto_genes_results,
          "@sample_name.@replicate_name.kallisto.isoforms.results" => :kallisto_isoforms_results,
          "@sample_name.@replicate_name.kallisto.genome.bam" => :kallisto_genome_bam,
          "@sample_name.@replicate_name.bwa_aligned.bam" => :bwa_aligned_bam,
          "@sample_name.@replicate_name.total.sam" => :total_sam,
          "@sample_name.@replicate_name.unique.sam" => :unique_sam
        },
        "@cohort_name" => {
          "@cohort_name.ug.raw.vcf" => :ug_raw_vcf,
          "@cohort_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@cohort_name.ug.filtered.vcf" => :ug_filtered_vcf,
          "@cohort_name.fpkm_table" => :fpkm_table,
          "@cohort_name.tpm_table" => :tpm_table,
          "@cohort_name.coverage_table" => :coverage_table,
          "@cohort_name.kallisto.coverage_table" => :kallisto_coverage_table
        }
      }
    })

    # tophat_align
    def_var :input_fastq1s do replicate.inputs.map(&:fq1) end
    def_var :input_fastq2s do replicate.inputs.map(&:fq2) end

    # count_transcripts
    def_var :transcripts_gtfs do replicates.map{ |r| transcripts_gtf r } end

    # assemble_transcripts
    def_var :gene_trackings do replicates.map{|r| gene_tracking(r) } end
    def_var :rsem_genes_resultses do replicates.map{|r| rsem_genes_results(r) } end
    def_var :transcripts_covs do replicates.map{|r| transcripts_cov(r) } end

    # qc
    def_var :qc_bam do replicate_bam end
    
    # bwa_qc
    def_var :bwa_qc_bam do replicate_bam end


    # diff_exp
    def_var :q_value_cutoff do 0.001 end
    
    def_var :mutations_config do "#{config_dir}/rna_mutations.yml" end

    def_var :find_novel_transcripts do nil end
  end
end
