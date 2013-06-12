module Exome
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end
    def platform_unit; "Exome"; end

    def_var :unit do "exome" end
    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do samples.map{ |s| sample_bam(s) } end
    def_var :bam_label do "bwa.realigned.dedup.recal" end
    def_var :input_bam do |s| s ? s.input_bam : nil end
    def_var :output_bams do samples.map{|s| output_bam(s) } end
    def_var :normal_bam do sample_bam(normal) end
    def_var :normal_name do sample.normal_name || samples.first.sample_name end
    def_var :normal do samples.find{|s| s.sample_name == normal_name} end
    def_var :tumor_bam do sample_bam(sample) end

    empty_var :cosmic_vcf

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "input.@input_name.read1.sai" => :read1_sai,
          "input.@input_name.read2.sai" => :read2_sai,
          "input.@input_name.mate_fixed.bwa.sam" => :mated_sam,
          "input.@input_name.paired.bwa.sam" => :paired_sam,
          "input.@input_name.renamed.bwa.sam" => :renamed_sam,
          "input.@input_name.aligned.bwa.bam" => :aligned_bam,

          "raw_sample.bam" => :raw_sample_bam,
          "raw_sample.bai" => :raw_sample_bai,

          "@chrom_name.snvs.raw.mutect.txt" => :mutect_snvs,
          "@chrom_name.snvs.raw.mutect.tmp.txt" => :mutect_snvs_tmp,
          "@sample_name.snvs.raw.mutect.txt" => :mutect_all_snvs,
          "@chrom_name.snvs.coverage.mutect.wig" => :mutect_coverage,
          "@chrom_name.insert_mutations" => :insert_mutations,
          "@chrom_name.snvs.pindel" => :pindel_snvs,
          "@chrom_name.snvs.pindel_D" => :pindel_snv_d,
          "@chrom_name.pindel.conf" => :pindel_list,
          "@chrom_name.indels.unpatched.pindel.vcf" => :pindel_unpatched_vcf,
          "@chrom_name.indels.raw.pindel.vcf" => :pindel_vcf,
          "@sample_name.indels.raw2.pindel.vcf" => :pindel_all_vcf,

          "@sample_name.ug.raw.vcf" => :ug_raw_vcf,
          "@sample_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@sample_name.ug.filtered.vcf" => :ug_filtered_vcf,

          ":normal_name.cov" => :normal_cov,
          "@sample_name.cov" => :tumor_cov,
          "absolute" => {
            "." => :absolute_scratch,
            "@sample_name.ABSOLUTE.RData" => :absolute_rdata
          }
        },
        "@cohort_name" => {
          "merged_library.bam" => :merged_library_bam,
          "raw_library.bam" => :raw_library_bam,
          "@chrom_name.merged.intervals" => :merged_intervals,
          "@chrom_name.realigned.bam" => :realigned_bam,
          "@chrom_name.recal.grp" => :recal_grp,
          "@chrom_name.recal.bam" => :recal_bam,
          "@chrom_name.recal.bai" => :recal_bai,
          "_splitbam_" => :split_bam_root,
          "_splitbam_@sample_name.bam" => :split_bam,
          "@cohort_name.interval_bed" => :interval_bed,
          "absolute" => {
            "." => :absolute_review_dir,
            "@cohort_name.PP-calls_tab.txt" => :review_table,
            "@cohort_name.PP-modes.RData" => :absolute_modes,
            "@cohort_name.PP-called_tab.txt" => :reviewed_table
          }
        }
      },
      ":metrics_dir" => {
        "@sample_name" => {
          "@sample_name.flagstat" => :qc_flag,
          "@sample_name.histogram" => :qc_histogram,
          "@sample_name.hybrid_selection_metrics" => :qc_hybrid,
          "@sample_name.alignment_metrics" => :qc_align_metrics,
          "@sample_name.insert_sizes" => :qc_inserts,
          "@sample_name.sample_summary" => :qc_coverage_metrics,
          "@sample_name" => :qc_coverage_base
        },
        "@cohort_name.duplication_metrics" => :duplication_metrics,
        "@cohort_name.qc_summary" => :qc_summary
      },
      ":output_dir" => {
        "@sample_name" => {
          "@sample_name.:bam_label.bam" => :output_bam,
          "@sample_name.mut.txt" => :tumor_muts,
          "@sample_name.mutations" => :sample_mutations,
          "@sample_name.gene_cnr" => :tumor_gene_cnr,
          "@sample_name.exon_cnr" => :tumor_exon_cnr,
          "@sample_name.cnr.Rdata" => :tumor_cnr_rdata,
          "@sample_name.cnr.seg" => :tumor_cnr_seg,
          "@sample_name.mutations" => :tumor_mutations,
          "@sample_name.normal_mut.txt" => :normal_muts,
        }
      }
    })

    def init_hook
      # add various bells and whistles here
      samples.each do |s|
        if s.inputs
          s.inputs.each do |i|
            i.add_member :input_name, i.index
          end
        end
        s.extend_with :chroms => chromosomes
      end
      @config.extend_with :chroms => chromosomes
    end

    def_var :chrom do job_item.chrom_name end

    # Align
    def_var :input_fastq1 do job_item.fq1 end
    def_var :input_fastq2 do job_item.fq2 end

    # Merge 
    def_var :aligned_bams do sample.inputs.map { |input| aligned_bam input } end
          
    # Library Merge
    def_var :raw_sample_bams do samples.map{|s| raw_sample_bam(s) } end

    # Library Split
    def_var :recal_bams do @config.chroms.map{ |c| recal_bam c } end
    def_var :split_bams do samples.map{ |s| split_bam s } end

    # Hybrid qc
    def_var :qc_bam do tumor_bam end

    #mut_filter
    def_var :pindel_vcfs do sample.chroms.map{|c| pindel_vcf c } end
    def_var :mutect_snvses do sample.chroms.map{|c| mutect_snvs c } end

    # Absolute
    def_var :absolute_rdatas do tumor_samples.map{|s| absolute_rdata s } end

    def_var :mutations_config do "#{config_dir}/exome_mutations.yml" end
  end
end
