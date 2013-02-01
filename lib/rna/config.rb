module Rna
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig
    include Pipeline::SampleConfig

    def_var :replicates do key_list :replicates end
    def_var :replicate_name do |i| "r#{i || replicate_index}" end
    def_var :replicate_scratch do "#{sample_scratch}/#{replicate_name}" end
    def_var :sample_names do key_list :sample_name end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do sample_names.map{|s| sample_bam(s) } end
    def_var :bam_label do "aligned.merged.sorted" end
    def_var :input_bam do |s| s ? sample(s)[:input_bam] : nil end
    def_var :output_bam do |s,r| "#{sample_output(s)}/#{s || sample_name}.#{r || replicate_name}.#{bam_label}.bam" end
    def_var :output_bams do samples.map{|s| s[:replicates].map_index{ |r,i| output_bam s[:sample_name], replicate_name(i) } }.flatten end
    def_var :normal_bam do sample_bam(normal_name) end
    def_var :normal_name do sample_names[0] end
    def_var :tumor_bam do sample_bam(sample_name) end

    # tophat_align
    def_var :tophat_scratch do "#{scratch}/tophat" end
    def_var :input_fastq1s do replicates[job_index][:input_fastq1_list] end
    def_var :input_fastq2s do replicates[job_index][:input_fastq2_list] end
    def_var :accepted_bam do "#{tophat_scratch}/accepted_hits.bam" end
    def_var :unmapped_bam do "#{tophat_scratch}/unmapped.bam" end
    def_var :merged_bam do "#{tophat_scratch}/merged.bam" end
    def_var :sorted_bam do "#{tophat_scratch}/merged.sorted.bam" end
    def_var :sorted_header do "#{tophat_scratch}/sorted_header.txt" end

    # count_transcripts
    def_var :cufflinks_scratch do "#{scratch}/cufflinks" end
    def_var :transcripts_gtf do "#{cufflinks_scratch}/transcripts.gtf" end
    def_var :output_gtf do "#{sample_output}/#{sample_name}.#{replicate_name}.transcripts.gtf" end
    def_var :transcripts_gtfs do samples.map{ |s| s[:replicates].map_index{ |r,i| "#{scratch_dir}/#{s[:sample_name]}/#{replicate_name(i)}/cufflinks/transcripts.gtf" } }.flatten end

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
