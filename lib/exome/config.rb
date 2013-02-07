module Exome
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig
    include Pipeline::SampleConfig

    def platform; "Illumina"; end
    def platform_unit; "Exome"; end

    def_var :unit do "exome" end
    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do sample_names.map{ |s| sample_bam(s) } end
    def_var :bam_label do "bwa.realigned.rmDups.recal" end
    def_var :input_bam do |s| s ? sample(s)[:input_bam] : nil end
    def_var :output_bam do |s| "#{output_dir}/#{s || sample_name}/#{s || sample_name}.#{bam_label}.bam" end
    def_var :output_bams do sample_names.map{|s| output_bam(s) } end
    def_var :normal_bam do sample_bam(normal_name) end
    def_var :normal_name do sample_names.first end
    def_var :tumor_bam do sample_bam(sample_name) end

    # Align
    def_var :input_fastq1 do input_fastq1_list.flatten[job_index] end
    def_var :input_fastq1_list do samples.map{ |s| s[:input_fastq1_list] }.flatten end
    def_var :input_fastq2 do input_fastq2_list.flatten[job_index] end
    def_var :input_fastq2_list do samples.map{ |s| s[:input_fastq2_list] }.flatten end

    def_var :sample_tag do key_index_list(:input_fastq1_list)[job_index] end
    def_var :scratch_tag_file do |affix,s,t| "#{scratch}/#{s || sample_name}.#{t || sample_tag}.#{affix}" end
    def_var :mated_sam do |s,t| scratch_tag_file "mateFixed.bwa.sam", s, t end
    def_var :read1_sai do |s,t| scratch_tag_file "read1.sai", s, t end
    def_var :read2_sai do |s,t| scratch_tag_file "read2.sai", s, t end
    def_var :paired_sam do |s,t| scratch_tag_file "paired.bwa.sam", s, t end
    def_var :renamed_sam do |s,t| scratch_tag_file "renamed.bwa.sam", s, t end
    def_var :renamed_bam do |s,t| scratch_tag_file "renamed.bwa.bam", s, t end


        # Merge stuff
    def_var :aligned_bams do samples[job_index][:input_fastq1_list].size.times.map { |t| renamed_bam sample_name, t } end
    def_var :merged_bam do |s| "#{sample_scratch(s)}/merged.bam" end
    def_var :merged_bai do "#{scratch}/merged.bam.bai" end
          
        # Recal stuff
    def_var :dedup_bam do "#{scratch}/dedup.bam" end
    def_var :mated_bam do "#{scratch}/mated.bam" end
    def_var :merged_intervals do "#{scratch}/merged.intervals" end
    def_var :realigned_bam do "#{scratch}/realigned.bam" end
    def_var :recal_bam do "#{scratch}/recal.bam" end
    def_var :recal_bai do "#{scratch}/recal.bam.bai" end
    def_var :recal_bams do sample_names.map{|s| merged_bam(s) } end
    def_var :remerged_bam do "#{scratch}/merged.bam" end
    def_var :recal_grp do "#{scratch}/recal.grp" end
    def_var :recal_metrics do "#{metrics_dir}/#{job_name}.recal_metrics" end
    def_var :split_bam do |s| s ? "#{scratch}/_splitbam_#{s}.bam" : "#{scratch}/_splitbam_" end
    def_var :split_bams do sample_names.map{ |s| split_bam(s) } end

        # hybrid_qc
    def_var :qc_bam do tumor_bam end
    def_var :qc_flag do |s| "#{metrics_dir}/#{s || sample_name}.flagstat" end
    def_var :qc_histogram do |s| "#{metrics_dir}/#{s || sample_name}.histogram" end
    def_var :qc_hybrid do |s| "#{metrics_dir}/#{s || sample_name}.hybrid_selection_metrics" end
    def_var :qc_align_metrics do |s| "#{metrics_dir}/#{s || sample_name}.alignment_metrics" end
    def_var :qc_inserts do |s| "#{metrics_dir}/#{s || sample_name}.insert_sizes" end

        # mut_det
    def_var :mutect_snvs do "#{output_dir}/#{sample_name}/#{sample_name}.snvs.raw.mutect.txt" end
    def_var :mutect_indels_anno do "#{output_dir}/#{sample_name}/#{sample_name}.indels.annotated.vcf" end
    def_var :mutect_indels_raw do "#{output_dir}/#{sample_name}/#{sample_name}.indels.raw.vcf" end
    def_var :mutect_indels_temp do "#{scratch}/#{sample_name}.indels.temp.vcf" end
    def_var :mutect_mutations do "#{scratch}/#{sample_name}.mutations" end
    def_var :mutect_coverage do "#{output_dir}/#{sample_name}/#{sample_name}.snvs.coverage.mutect.wig" end
    def_var :insert_mutations do "#{scratch}/#{sample_name}.insert_mutations" end


        # copy_number
    def_var :interval_bed do "#{job_scratch}/intervals.bed" end
    def_var :normal_cov do "#{scratch}/#{normal_name}.cov" end
    def_var :tumor_cov do "#{scratch}/#{sample_name}.cov" end
    def_var :tumor_gene_cnr do "#{output_dir}/#{sample_name}/#{sample_name}.gene_cnr" end
    def_var :tumor_exon_cnr do "#{output_dir}/#{sample_name}/#{sample_name}.exon_cnr" end
    def_var :tumor_cnr_rdata do "#{output_dir}/#{sample_name}/#{sample_name}.cnr.Rdata" end
    def_var :tumor_cnr_seg do "#{output_dir}/#{sample_name}/#{sample_name}.cnr.seg" end
    def_var :tumor_mutations do "#{output_dir}/#{sample_name}/#{sample_name}.mutations" end
  end
end
