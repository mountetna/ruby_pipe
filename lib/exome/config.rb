module Exome
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end
    def platform_unit; "Exome"; end

    def_var :unit do "exome" end
    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do sample_names.map{ |s| sample_bam(s) } end
    def_var :bam_label do "bwa.realigned.dedup.recal" end
    def_var :input_bam do |s| s ? sample(s)[:input_bam] : nil end
    def_var :output_bam do |s| sample_output_file "#{s || sample_name}.#{bam_label}.bam", s end
    def_var :output_bams do sample_names.map{|s| output_bam(s) } end
    def_var :normal_bam do sample_bam(normal_name) end
    def_var :normal_name do sample.normal_name || sample_names.first end
    def_var :tumor_bam do sample_bam(sample_name) end

    # Align
    def_var :input_fastq1 do job_item.fq1 end
    def_var :input_fastq2 do job_item.fq2 end
    def_var :align_tag do |i| (i || job_item).index end
    def_var :align_file do |affix,s,i| sample_scratch_file "input.#{align_tag(i)}.#{affix}", s end
    def_var :read1_sai do align_file "read1.sai" end
    def_var :read2_sai do align_file "read2.sai" end
    def_var :mated_sam do align_file "mate_fixed.bwa.sam" end
    def_var :paired_sam do align_file "paired.bwa.sam" end
    def_var :renamed_sam do align_file "renamed.bwa.sam" end
    def_var :aligned_bam do |s,i| align_file "aligned.bwa.bam", s, i end

    # Merge 
    def_var :aligned_bams do sample.inputs.map { |input| aligned_bam sample_name, input } end
    def_var :raw_sample_bam do |s| sample_scratch_file "raw_sample.bam", s end
    def_var :raw_sample_bai do sample_scratch_file "raw_sample.bam.bai" end
          
    # Library Merge
    def_var :raw_sample_bams do sample_names.map{|s| raw_sample_bam(s) } end
    def_var :merged_library_bam do cohort_scratch_file "merged_library.bam" end
    def_var :raw_library_bam do cohort_scratch_file "raw_library.bam" end
    def_var :recal_metrics do cohort_metrics_file "recal.metrics" end

    # Recal
    def_var :chrom do job_item[:contig] end
    def_var :chrom_file do |affix,c| cohort_scratch_file "#{c || chrom}.#{affix}" end
    def_var :merged_intervals do chrom_file "merged.intervals" end
    def_var :realigned_bam do chrom_file "realigned.bam" end
    def_var :recal_grp do chrom_file "recal.grp" end
    def_var :recal_bam do |c| chrom_file "recal.bam", c end
    def_var :recal_bai do chrom_file "recal.bam.bai" end

    # Library Split
    def_var :recal_bams do chroms.map{ |c| recal_bam c[:contig] } end
    def_var :split_bam_root do |e| cohort_scratch_file "_splitbam_#{e}" end
    def_var :split_bam do |s| split_bam_root(s || sample_name) + ".bam" end
    def_var :split_bams do sample_names.map{ |s| split_bam(s) } end

    # Hybrid qc
    def_var :qc_bam do tumor_bam end
    def_var :qc_flag do |s| sample_metrics_file "flagstat", s end
    def_var :qc_histogram do |s| sample_metrics_file "histogram", s end
    def_var :qc_hybrid do |s| sample_metrics_file "hybrid_selection_metrics", s end
    def_var :qc_align_metrics do |s| sample_metrics_file "alignment_metrics", s end
    def_var :qc_inserts do |s| sample_metrics_file "insert_sizes", s end
    def_var :qc_coverage_metrics do |s| sample_metrics_file "sample_summary", s end

    # mut_det
    def_var :mutect_snvs do |c| sample_scratch_file "#{c||chrom}.snvs.raw.mutect.txt" end
    def_var :sample_mutations do sample_output_file "mutations" end
    def_var :mutect_coverage do sample_scratch_file "#{chrom}.snvs.coverage.mutect.wig" end
    def_var :insert_mutations do sample_scratch_file "#{chrom}.insert_mutations" end
    def_var :pindel_snvs do sample_scratch_file "#{chrom}.snvs.pindel" end
    def_var :pindel_snv_d do sample_scratch_file "#{chrom}.snvs.pindel_D" end
    def_var :pindel_list do sample_scratch_file "#{chrom}.pindel.conf" end
    def_var :pindel_vcf do |c| sample_scratch_file "#{c || chrom}.indels.raw.pindel.vcf" end

    #mut_filter
    def_var :pindel_vcfs do chroms.map{|c| pindel_vcf c[:contig] } end
    def_var :mutect_snvses do chroms.map{|c| mutect_snvs c[:contig] } end

    # copy_number
    def_var :interval_bed do cohort_scratch_file "intervals.bed" end
    def_var :normal_cov do sample_scratch_file "#{normal_name}.cov" end
    def_var :tumor_cov do sample_scratch_file "#{sample_name}.cov" end
    def_var :tumor_gene_cnr do sample_output_file "gene_cnr" end
    def_var :tumor_exon_cnr do sample_output_file "exon_cnr" end
    def_var :tumor_cnr_rdata do sample_output_file "cnr.Rdata" end
    def_var :tumor_cnr_seg do sample_output_file "cnr.seg" end
    def_var :tumor_mutations do sample_output_file "mutations" end

    def_var :mutations_config do "#{config_dir}/exome_mutations.yml" end
  end
end
