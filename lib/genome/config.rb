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
          "input.@input_name.chaste.bam" => :chaste_bam,
          "input.@input_name.reads1.fastq.gz" => :reads1_fastq,
          "input.@input_name.reads2.fastq.gz" => :reads2_fastq,

          "chunk_info" => :chunk_info,

          "input.read1.sai" => :read1_sai,
          "input.read2.sai" => :read2_sai,
          "@input_name.paired.bwa.sam" => :paired_sam,

          "chunk.@chunk_name.read1.fastq.gz" => :chunk1_fastq,
          "chunk.@chunk_name.read2.fastq.gz" => :chunk2_fastq,
          "chunk.@chunk_name.read1.sai" => :read1_sai,
          "chunk.@chunk_name.read2.sai" => :read2_sai,
          "chunk.@chunk_name.mate_fixed.bwa.bam" => :mated_bam,
          "chunk.@chunk_name.paired.bwa.sam" => :paired_sam,
          "chunk.@chunk_name.dedup.bwa.bam" => :dedup_bam,
          "chunk.@chunk_name.aligned.bwa.bam" => :aligned_bam,

          "raw_sample.bam" => :raw_sample_bam,
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
          "@sample_name.recal.bam" => :recal_bam,
          "absolute" => {
            "." => :absolute_scratch,
            "@sample_name.ABSOLUTE.RData" => :absolute_rdata
          }
        },
        "@lane_name" => {
          "merged_lane.bam" => :merged_lane_bam,
          "recal.grp" => :recal_grp,
          "recal_plot.pdf" => :recal_plot_pdf,
          "_splitbam_" => :lane_split_bam_root,
          "_splitbam_SAMPLE.bam" => :lane_split_bam
        },
        "@patient_name" => {
          "merged_patient.bam" => :merged_patient_bam,
          "raw_patient.bam" => :raw_patient_bam,
          "duplication_metrics" => :duplication_metrics,
          "@chrom_name.patient.intervals" => :patient_intervals,
          "@chrom_name.realigned_patient.bam" => :realigned_patient_bam,
          "@chrom_name.split." => :patient_split_bam_root,
          "@chrom_name.split.@sample_name.bam" => :patient_split_bam,
          "@patient_name.ug.raw.vcf" => :ug_raw_vcf,
          "@patient_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@patient_name.ug.filtered.vcf" => :ug_filtered_vcf
        },
        "@cohort_name" => {
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

    def make_chunks s
      if File.exists? chunk_info(s)
        num_reads = File.read(chunk_info(s)).to_i
        num_chunks = (num_reads*4 / chunk_size.to_f).ceil
        num_chunks.times.map{ |i|
          { :chunk_number => i, :chunk_name => "c#{i}" }
        }
      else
        # you don't know right now, just make a dummy
        [ { :chunk_number => 0, :chunk_name => "c0" } ]
      end
    end

    def make_lanes
      samples.group_by(&:lane_name).map do |lane,samples|
        { :lane_name => lane, :samples => samples }
      end
    end

    def make_patients
      samples.group_by(&:patient_name).map do |patient,samples|
        { :patient_name => patient, :samples => samples }
      end
    end

    def make_fastq filename
      {
        :fastq_name => File.basename(filename).sub(/.fastq.gz/,""),
        :fastq_file => filename
      }
    end

    def init_hook
      # add various bells and whistles here
      samples.each do |s|
        if s.inputs
          s.inputs.each do |i|
            i.add_member :input_name, i.index
          end
          s.extend_with :fastqs => s.inputs.map{|i|
            [ make_fastq(i.fq1), make_fastq(i.fq2) ]
          }.flatten
        end
        s.extend_with :chroms => chromosomes
        s.extend_with :chunks => make_chunks(s)
        s.add_member :lane_name, "lane#{s.lane || 0}"
        s.add_member :patient_name, "patient#{s.patient || 0}"
      end
      @config.extend_with :lanes => make_lanes
      lanes.each do |lane|
        lane.extend_with :chroms => chromosomes
      end
      @config.extend_with :patients => make_patients
      patients.each do |patient|
        patient.extend_with :chroms => chromosomes
      end
      @config.extend_with :chroms => chromosomes
    end

    def_var :initial_bam do |s| (s || job_item).initial_bam end
    
    def_var :lane_name do |l| (l||job_item).property :lane_name end
    def_var :chunk_size do 4_000_000 end
    job_items :chrom, :lane, :patient, :chunk, :fastq

    # Align
    def_var :reads_bam do job_item.reads_bam end
    def_var :paired_sams do sample.inputs.map{|input| paired_sam input } end
    def_var :input_fastq1 do |i| (i || job_item).fq1 end
    def_var :input_fastq2 do |i| (i || job_item).fq2 end

    def_var :reads1_fastqs do |s| (s || sample).inputs.map{|input| input_fastq1(input) || reads1_fastq(input) } end
    def_var :reads2_fastqs do |s| (s || sample).inputs.map{|input| input_fastq2(input) || reads2_fastq(input) } end
 
    # Merge 
    def_var :aligned_bams do sample.chunks.map { |chunk| aligned_bam chunk } end
          
    # Lane Merge
    def_var :lane_raw_sample_bams do lane.samples.map{|s| raw_sample_bam(s) } end

    # Patient Merge
    def_var :patient_raw_sample_bams do patient.samples.map {|s| find_prop :patient_raw_sample_bam, s } end
    def_var :patient_raw_sample_bam do |s| find_prop :recal_bam, s end

    # PatientSplit
    def_var :realigned_patient_bams do patient.chroms.map{|c| realigned_patient_bam(c) } end
    def_var :patient_split_bams do patient.samples.map{|s| patient_split_bam s.chroms.find{|c| c.chrom_name == chrom.chrom_name} } end
    def_var :patient_sample_bams do patient.samples.map{|s| sample_bam(s) } end
    def_var :sample_patient_bams do sample.chroms.map{|c| patient_split_bam c } end

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
