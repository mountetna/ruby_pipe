module Rna
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig
    include Pipeline::SampleConfig

    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do sample_names.map{|s| sample_bam(s) } end
    def_var :bam_label do "aligned.merged.sorted" end
    def_var :input_bam do |s| s ? sample(s)[:input_bam] : nil end
    def_var :output_bam do |s,r| sample_output_file  "#{s || sample_name}.#{r || replicate_name}.#{bam_label}.bam", s end
    def_var :output_bams do replicates.map{ |r| output_bam r.parent.sample_name, replicate_name(r.index) } end
    def_var :normal_bam do sample_bam(normal_name) end
    def_var :normal_name do sample.normal_name || sample_names.first end
    def_var :tumor_bam do sample_bam(sample_name) end

    def_var :replicate do job_item end
    def_var :replicate_index do job_index end
    def_var :replicates do samples.collect &:replicates end
    def_var :replicate_name do |i| "r#{i || replicate_index}" end
    def_var :replicate_scratch do |r| sample_scratch_dir(r || replicate_name) end
    def_var :replicate_scratch_file do |affix,r| File.join replicate_scratch(r), affix end

    # tophat_align
    def_var :tophat_scratch do replicate_scratch_file "tophat" end
    def_var :tophat_scratch_file do |affix| File.join tophat_scratch, affix end
    def_var :input_fastq1s do replicate.inputs.map(&:fq1) end
    def_var :input_fastq2s do replicate.inputs.map(&:fq2) end
    def_var :accepted_bam do tophat_scratch_file "accepted_hits.bam" end
    def_var :unmapped_bam do tophat_scratch_file "unmapped.bam" end
    def_var :merged_bam do tophat_scratch_file "merged.bam" end
    def_var :sorted_bam do tophat_scratch_file "merged.sorted.bam" end
    def_var :sorted_header do tophat_scratch_file "sorted_header.txt" end

    # count_transcripts
    def_var :cufflinks_scratch do |r| replicate_scratch_file "cufflinks" end
    def_var :cufflinks_scratch_file do |affix,r| File.join cufflinks_scratch, affix end
    def_var :transcripts_gtf do cufflinks_scratch_file "transcripts.gtf" end
    def_var :output_gtf do sample_output_file "#{sample_name}.#{replicate_name}.transcripts.gtf" end
    def_var :transcripts_gtfs do
      replicates.map{ |r| cufflinks_scratch_file
      samples.map{ |s| s[:replicates].map_index{ |r,i| "#{scratch_dir}/#{s[:sample_name]}/#{replicate_name(i)}/cufflinks/transcripts.gtf" } }.flatten end

    # assemble_transcripts
    def_var :cuffcompare_scratch do "#{scratch}/cuffcompare/cuffcomp" end
    def_var :tracking_file do "#{cuffcompare_scratch}.tracking" end
    def_var :assembly_list do "#{scratch}/assembly.list" end
    def_var :assembly_gtf do "#{scratch}/merged.gtf" end

    # qc
    def_var :qc_bam do tumor_bam end
    def_var :qc_flag do |s| "#{metrics_dir}/#{s || sample_name}.#{replicate_name}.flagstat" end
    def_var :qc_rnaseq do |s| "#{metrics_dir}/#{s || sample_name}.#{replicate_name}.rnaseq_metrics" end
    def_var :qc_pdf do |s| "#{metrics_dir}/#{s || sample_name}.#{replicate_name}.rnaseq.pdf" end
    def_var :qc_align_metrics do |s| "#{metrics_dir}/#{s || sample_name}.#{replicate_name}.alignment_metrics" end

    #univ_geno
    def_var :ug_raw_vcf do "#{output_dir}/#{job_name}/#{job_name}.ug.raw.vcf" end
    def_var :ug_annotated_vcf do "#{scratch}/#{job_name}.ug.annotated.vcf" end
    def_var :ug_filtered_vcf do "#{output_dir}/#{job_name}/#{job_name}.ug.filtered.vcf" end

    #filter_muts
    def_var :sample_mutations do "#{sample_output}/#{sample_name}.mutations" end
    def_var :mutations_config do "#{config_dir}/rna_mutations.yml" end
  end
end
