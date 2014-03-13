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

    def_var :bam_label do "aligned.sorted" end
    def_var :fdr_cutoff do 0.10 end
    def_var :normal_bams do replicate_bams normal end

    def_var :replicate_bams do |s| (s||sample).replicates.map{ |r| replicate_bam r} end
    def_var :replicate_bam do |r| input_bam(r) || output_bam(r) end
    def_var :replicates do |s| (s||sample).replicates end

    def_var :diff_exps do samples.collect(&:diff_exp).compact.flatten end

    def_var :sample_replicate_name do |r| "#{r.property :sample_name}.#{r.property :replicate_name}" end

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "@replicate_name" => {
            "coverage.sam" => :coverage_sam,
            "splice.sam" => :splice_sam,
            "chimera1.fastq" => :chimera_fq1,
            "chimera2.fastq" => :chimera_fq2,
            "chimerascan" => {
              "." => :chimera_dir,
              "chimeras.bedpe" => :chimera_bedpe
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
              "tmp" => :rsem_tmp_dir,
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
            }
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
          "@sample_name.@replicate_name.alignment_metrics" => :qc_align_metrics
        }
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
          "@sample_name.@replicate_name.exon_splice_counts" => :exon_splice_counts
        },
        "@cohort_name" => {
          "@cohort_name.ug.raw.vcf" => :ug_raw_vcf,
          "@cohort_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@cohort_name.ug.filtered.vcf" => :ug_filtered_vcf,
          "@cohort_name.fpkm_table" => :fpkm_table,
          "@cohort_name.tpm_table" => :tpm_table,
          "@cohort_name.coverage_table" => :coverage_table
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

    # diff_exp
    def_var :q_value_cutoff do 0.001 end
    
    def_var :mutations_config do "#{config_dir}/rna_mutations.yml" end

    def_var :find_novel_transcripts do nil end
  end
end
