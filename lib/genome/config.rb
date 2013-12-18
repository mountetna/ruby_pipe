module Genome
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end
    def platform_unit; "Genome"; end

    def_var :unit do "genome" end
    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do samples.map{ |s| sample_bam(s) } end
    def_var :bam_label do "bwa.realigned.dedup.recal" end
    def_var :input_bam do |s| (s || sample).input_bam end
    def_var :output_bams do samples.map{|s| output_bam(s) } end
    def_var :normal_bam do sample_bam(normal) end
    def_var :normal_name do sample.normal_name || samples.first.sample_name end
    def_var :normal do samples.find{|s| s.sample_name == normal_name} end
    def_var :tumor_bam do sample_bam(sample) end

    empty_var :cosmic_vcf

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "input.read1.sai" => :read1_sai,
          "input.read2.sai" => :read2_sai,
          "@input_name.paired.bwa.sam" => :paired_sam,

          "paired_merge.sam" => :paired_merge_sam,
          "paired_merge.bam" => :paired_merge_bam,
          "mate_fixed.bam" => :mated_bam,
          "raw_sample.bam" => :raw_sample_bam,
          "dedup_sample.bam" => :dedup_sample_bam,
          "raw_sample.bai" => :raw_sample_bai,

          "@chrom_name.normal.snp.txt" => :normal_mut,

          "@sample_name.snps.baf"=> :tumor_baf,
          ":normal_name.snps.baf"=> :normal_baf,
          "@chrom_name.snvs.raw.mutect.txt" => :mutect_snvs,
          "@sample_name.@chrom_name.snvs.joint.raw.variant.vcf" => :snp_vcf,
          "@sample_name.@chrom_name.snvs.snp.annotated.vcf" => :snp_annotated_vcf,
          "@sample_name.@chrom_name.snvs.snp.filtered.vcf" => :snp_filtered_vcf,
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
          "@sample_name.ratio" => :tumor_ratio,
          "@sample_name.all_muts.maf" => :all_muts_maf,
          "@sample_name.pindel.output.txt" => :output_for_pindel,
          "@sample_name.temp.pindel.output.txt" => :temp_output_for_pindel,

          ":normal_name.cov" => :normal_cov,
          "@sample_name.cov" => :tumor_cov,
          "absolute" => {
            "." => :absolute_scratch,
            "@sample_name.ABSOLUTE.RData" => :absolute_rdata
          }
        },
        "@cohort_name" => {
          #for copy number
          ":normal_name.:bam_label.@chrom_name" => :split_normal_contigs,
          "@sample_name.:bam_label.@chrom_name" => :split_tumor_contigs,

          "merged_library.bam" => :merged_library_bam,
          "raw_library.bam" => :raw_library_bam,
          "@chrom_name.merged.intervals" => :merged_intervals,
          "@chrom_name.realigned.bam" => :realigned_bam,
          "@chrom_name.recal.grp" => :recal_grp,
	  "@chrom_name.post.recal.grp" => :post_recal_grp,
	  "@chrom_name.pre.recal.pdf" => :pre_recal_plot,
	  "@chrom_name.post.recal.pdf" => :post_recal_plot,
          "@chrom_name.recal.bam" => :recal_bam,
          "@chrom_name.recal.bai" => :recal_bai,
          "_splitbam_" => :split_bam_root,
	  "@chrom_name.sample.bam" => :split_sample_root,
          "_splitbam_@sample_name.bam" => :split_bam,
          "_splitsamplebam_@sample_name.bam" => :split_sample_bam,
          "@cohort_name.bed" => :interval_bed,
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
          "@sample_name.alignment_metrics" => :qc_align_metrics,
          "@sample_name.insert_sizes" => :qc_inserts,
          "@sample_name.sample_summary" => :qc_coverage_metrics,
          "@sample_name" => :qc_coverage_base,
          "@sample_name.duplication_metrics" => :duplication_metrics,
        },
        "@cohort_name.qc_summary" => :qc_summary
      },
      ":output_dir" => {
        "@sample_name" => {
          "@sample_name.:bam_label.bam" => :output_bam,
          "@sample_name.mut.txt" => :tumor_muts,
          "@sample_name.somatic.maf" => :tumor_maf,
          "@sample_name.germline.maf" => :germline_maf,
          "@sample_name.mutations" => :sample_mutations,
          "@sample_name.snp.filtered.annotated.vcf" => :ug_filtered_vcf,
          "@sample_name.gene_cnr" => :tumor_gene_cnr,
          "@sample_name.exon_cnr" => :tumor_exon_cnr,
          "@{sample_name}__.Rdata" => :tumor_cnr_rdata,
          "@{sample_name}_pscbs.Rdata" => :tumor_cnr_rdata_pscbs,
          "@sample_name.pscbs.cnr.seg" => :tumor_cnr_seg_pscbs,
          "@sample_name.cnr.seg" => :tumor_cnr_seg,
          "@sample_name.mutations" => :tumor_mutations,
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

    def_var :initial_bam do |s| (s || job_item).initial_bam end
    
    def_var :chrom do job_item.chrom_name end

    # Align
    def_var :input_fastq1 do job_item.fq1 end
    def_var :input_fastq2 do job_item.fq2 end
    def_var :paired_sams do sample.inputs.map{|input| paired_sam input } end
 
    # Merge 
    def_var :aligned_bams do sample.inputs.map { |input| aligned_bam input } end
          
    # Library Merge
    def_var :raw_sample_bams do samples.map{|s| raw_sample_bam(s) } end
    def_var :dedup_sample_bams do samples.map{|s| dedup_sample_bam(s) } end

    # Library Split
    def_var :recal_bams do @config.chroms.map{ |c| recal_bam c } end
    def_var :split_bams do samples.map{ |s| split_bam s } end

    #copy number
    def_var :tumor_cnr_rdatas do sample.chroms.map{ |c| tumor_cnr_rdata c } end
    
    # Hybrid qc
    def_var :qc_bam do tumor_bam end

    # Absolute
    def_var :absolute_rdatas do tumor_samples.map{|s| absolute_rdata s } end

    #mut_filter
    def_var :output_for_pindels do samples.map{ |s| output_for_pindel(s)} end
    def_var :temp_output_for_pindels do samples.map{ |s| temp_output_for_pindel s} end
    def_var :pindel_vcfs do sample.chroms.map{|c| pindel_vcf c } end
    def_var :mutect_snvses do sample.chroms.map{|c| mutect_snvs c } end
    def_var :snp_filtered_vcfs do sample.chroms.map{|c| snp_filtered_vcf c } end
    def_var :mutations_config do "#{config_dir}/genome_mutations.yml" end
  end
end
