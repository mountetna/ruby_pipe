module Exome
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end
    def platform_unit; "Exome"; end

    def_var :unit do "exome" end
    def_var :bam_label do "bwa.realigned.dedup.recal" end
    def_var :output_bams do samples.map{|s| output_bam(s) } end

    empty_var :cosmic_vcf

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "input.@input_name.chaste.bam" => :chaste_bam,
          "input.@input_name.reads1.fastq.gz" => :reads1_fastq,
          "input.@input_name.reads2.fastq.gz" => :reads2_fastq,

          "chunk.@chunk_name.read1.fastq.gz" => :chunk1_fastq,
          "chunk.@chunk_name.read2.fastq.gz" => :chunk2_fastq,
          "chunk.@chunk_name.read1.sai" => :read1_sai,
          "chunk.@chunk_name.read2.sai" => :read2_sai,
          "chunk.@chunk_name.mate_fixed.bwa.bam" => :mated_bam,
          "chunk.@chunk_name.paired.bwa.sam" => :paired_sam,
          "chunk.@chunk_name.paired.bwa.bam" => :paired_bam,
          "chunk.@chunk_name.aligned.bwa.bam" => :aligned_bam,

          "@sample_name.reads1.fastq" => :sample_reads1_fastq,
          "@sample_name.reads2.fastq" => :sample_reads2_fastq,

          "chunk_info" => :chunk_info,

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
          "@sample_name.recal.bam" => :recal_sample_bam,
          "absolute" => {
            "." => :absolute_scratch,
            "@sample_name.ABSOLUTE.RData" => :absolute_rdata
          }
        },
        "@lane_name" => {
          "merged_lane.bam" => :merged_lane_bam,
          "recal.grp" => :recal_grp,
          "recal_plot.pdf" => :recal_plot_pdf,
          "@chrom_name.recal.bam" => :recal_bam,
          "_splitbam_" => :lane_split_bam_root,
          "_splitbam_SAMPLE.bam" => :lane_split_bam
        },
        "@patient_name" => {
          "merged_patient.bam" => :merged_patient_bam,
          "raw_patient.bam" => :raw_patient_bam,
          "duplication_metrics" => :duplication_metrics,
          "@chrom_name.patient.intervals" => :patient_intervals,
          "@chrom_name.realigned_patient.bam" => :realigned_patient_bam,
          "_splitbam_" => :patient_split_bam_root,
          "_splitbam_SAMPLE.bam" => :patient_split_bam
        },
        "@cohort_name" => {
          "@chrom_name.realigned.bam" => :realigned_bam,
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
          "@sample_name.hybrid_selection_metrics" => :qc_hybrid,
          "@sample_name.alignment_metrics" => :qc_align_metrics,
          "@sample_name.insert_sizes" => :qc_inserts,
          "@sample_name.sample_summary" => :qc_coverage_metrics,
          "@sample_name.indelocator" => :indelocator_metrics,
          "@sample_name" => :qc_coverage_base
        },
        "@cohort_name.duplication_metrics" => :duplication_metrics,
        "@cohort_name.qc_summary" => :qc_summary
      },
      ":output_dir" => {
        "@sample_name" => {
          "@sample_name.:bam_label.bam" => :output_bam,
          "@sample_name.:bam_label.bam.bai" => :sample_bai,
          "@sample_name.mut.txt" => :tumor_muts,
          "@sample_name.maf" => :tumor_maf,
          "@sample_name.mutations" => :sample_mutations,
          "@sample_name.gene_cnr" => :tumor_gene_cnr,
          "@sample_name.exon_cnr" => :tumor_exon_cnr,
          "@sample_name.cnr.Rdata" => :tumor_cnr_rdata,
          "@sample_name.cnr.seg" => :tumor_cnr_seg,
          "@sample_name.mutations" => :tumor_mutations,
          "@sample_name.normal_mut.txt" => :normal_muts,
          "@sample_name.indelocator.bed" => :indelocator_bed,
          "@sample_name.indelocator.txt" => :indelocator_output,
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

    def init_hook
      # add various bells and whistles here
      samples.each do |s|
        if s.inputs
          s.inputs.each do |i|
            i.add_member :input_name, i.index
          end
        end
        s.extend_with :chroms => chromosomes
        s.extend_with :chunks => make_chunks(s)
      end
      @config.extend_with :lanes => make_lanes
      lanes.each do |lane|
        lane.extend_with :chroms => chromosomes
      end
      @config.extend_with :patients => make_patients
      patients.each do |patient|
        patient.extend_with :chroms => chromosomes
      end
    end

    def make_lanes
      samples.group_by(&:lane).map do |lane,samples|
        { :lane_name => "lane#{lane || 0}", :samples => samples }
      end
    end

    def make_patients
      samples.group_by(&:patient).map do |patient,samples|
        { :patient_name => "patient#{patient || 0}", :samples => samples }
      end
    end

    def_var :lane_name do |l| (l||job_item).property :lane_name end
    def_var :chunk_size do 4_000_000 end
    job_items :chrom, :lane, :patient

    def_var :reads1_fastqs do sample.inputs.map{|input| input_fastq1(input) || reads1_fastq(input) } end
    def_var :reads2_fastqs do sample.inputs.map{|input| input_fastq2(input) || reads2_fastq(input) } end

    # Align
    def_var :input_file1 do input_fastq1 || reads_bam end
    def_var :input_file2 do input_fastq2 || reads_bam end
    def_var :reads_bam do job_item.reads_bam end
    def_var :input_fastq1 do |i| (i || job_item).fq1 end
    def_var :input_fastq2 do |i| (i || job_item).fq2 end

    # Merge 
    def_var :aligned_bams do sample.chunks.map { |chunk| aligned_bam chunk } end
          
    # Lane Merge
    def_var :lane_raw_sample_bams do lane.samples.map{|s| raw_sample_bam(s) } end

    # Lane Split
    def_var :lane_recal_bams do lane.chroms.map{|c| recal_bam(c) } end
    def_var :lane_sample_split_bam do |s| lane_split_bam.sub("SAMPLE",s.sample_name) end
    def_var :lane_split_bams do lane.samples.map{|s| lane_sample_split_bam(s) } end
    def_var :lane_sample_bams do lane.samples.map{|s| recal_sample_bam(s) } end

    # Patient Merge
    def_var :patient_raw_sample_bams do patient.samples.map {|s| patient_raw_sample_bam s } end
    def_var :patient_raw_sample_bam do |s| recal_sample_bam(s || job_item) end

    # PatientSplit
    def_var :realigned_patient_bams do patient.chroms.map{|c| realigned_patient_bam(c) } end
    def_var :patient_sample_split_bam do |s| patient_split_bam.sub("SAMPLE",s.sample_name) end
    def_var :patient_split_bams do patient.samples.map{|s| patient_sample_split_bam(s) } end
    def_var :patient_sample_bams do patient.samples.map{|s| sample_bam(s) } end

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
