module Exome
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def platform; "Illumina"; end
    def platform_unit; "Exome"; end

    def_var :unit do "exome" end
    def_var :bam_label do "bwa.realigned.dedup.recal" end
    def_var :output_bams do samples.map{|s| output_bam(s) } end

    empty_var :cosmic_vcf, :verify_fastq_quality

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "@{fastq_name}_fastqc" => {
            "." => :fastqc_output_dir,
            "summary.txt" => :fastqc_summary,
            "fastqc_report.html" => :fastqc_html
          },

          "input.@input_name.chaste.bam" => :chaste_bam,
          "input.@input_name.chaste.sorted.bam" => :chaste_sorted_bam,
          "input.@input_name.reads1.fastq.gz" => :reads1_fastq,
          "input.@input_name.reads2.fastq.gz" => :reads2_fastq,

          "chunk.@chunk_name.read1.fastq.gz" => :chunk1_fastq,
          "chunk.@chunk_name.read2.fastq.gz" => :chunk2_fastq,
          "chunk.@chunk_name.read1.sai" => :read1_sai,
          "chunk.@chunk_name.read2.sai" => :read2_sai,
          "chunk.@chunk_name.mate_fixed.bwa.bam" => :mated_bam,
          "chunk.@chunk_name.paired.bwa.sam" => :paired_sam,
          "chunk.@chunk_name.dedup.bwa.bam" => :dedup_bam,
          "chunk.@chunk_name.aligned.bwa.bam" => :aligned_bam,

          "@sample_name.reads1.fastq" => :sample_reads1_fastq,
          "@sample_name.reads2.fastq" => :sample_reads2_fastq,

          "chunk_info" => :chunk_info,

          "raw_sample.bam" => :raw_sample_bam,
          "raw_sample.bai" => :raw_sample_bai,

          "@chrom_name.snvs.raw.mutect.txt" => :mutect_snvs,
          "@chrom_name.snvs.raw.mutect.vcf" => :mutect_vcf,
          "@chrom_name.snvs.annotated.mutect.vcf" => :mutect_annotated_vcf,
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
          "@chrom_name.indels.raw.strelka.vcf" => :strelka_vcf,
          "@chrom_name.strelka.realigned.bam" => :strelka_bam,
          "@chrom_name.indelocator.bed" => :indelocator_bed,
          "@chrom_name.indelocator.txt" => :indelocator_output,
          "@chrom_name.indelocator" => :indelocator_metrics,

          "@chrom_name.somaticindel.unpatched.vcf" => :somaticindel_unpatched_vcf,
          "@chrom_name.somaticindel.verbose.txt" => :somaticindel_verbose,
          "@chrom_name.somaticindel.indels.raw.vcf" => :somaticindel_vcf,
          "@chrom_name.somaticindel.indels.annotated.vcf" => :somaticindel_annotated_vcf,

          "@chrom_name.somatic.maf" => :tumor_chrom_maf,
          "@chrom_name.germline.maf" => :germline_chrom_maf,
          "@chrom_name.germline.unfiltered.vcf" => :germline_unfiltered_vcf,
          "@chrom_name.germline.annotated.vcf" => :germline_annotated_vcf,
          "@chrom_name.all_muts.maf" => :all_muts_chrom_maf,

          "@sample_name.context.maf" => :tumor_context_maf,

          "@sample_name.ug.SNPs.vcf" => :ug_snps_vcf,

          "@sample_name.tumor.baf" => :tumor_baf,
          "@sample_name.normal.baf" => :normal_baf,

          "@sample_name.cov" => :sample_cov,
          "@sample_name.cov_gc" => :sample_cov_gc,
          "@sample_name.cnvkit.target.cov" => :sample_target_cov,
          "@sample_name.cnvkit.reference.cnn" => :sample_reference_cnn,
          "@sample_name.cnvkit.antitarget.cov" => :sample_antitarget_cov,
          "@sample_name.recal.bam" => :recal_bam,
          "absolute" => {
            "." => :absolute_scratch,
            "@sample_name.ABSOLUTE.RData" => :absolute_scratch_rdata,
            "@sample_name.ABSOLUTE_plot.pdf" => :absolute_scratch_pdf
          }
        },
        ":normal_name" => {
          ":normal_name.cov" => :normal_cov,
          ":normal_name.cnvkit.reference.cnn" => :normal_reference_cnn,
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
          "_splitbam_@sample_name.bam" => :sample_split_bam,
          "@patient_name.ug.raw.vcf" => :ug_raw_vcf,
          "@patient_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@patient_name.ug.filtered.vcf" => :ug_filtered_vcf
        },
        "@cohort_name" => {
          "@chrom_name.realigned.bam" => :realigned_bam,
          "@cohort_name.bed" => :interval_bed,
          "@cohort_name.total_normal.cov" => :total_normal_cov,
          "absolute" => {
            "." => :absolute_review_dir,
            "@cohort_name.PP-calls_tab.txt" => :review_table,
            "@cohort_name.PP-modes.data.RData" => :absolute_modes,
            "@cohort_name.PP-called_tab.txt" => :reviewed_table,
            "reviewed" => {
              "@cohort_name.Exome.Pipeline.ABSOLUTE.table.txt" => :absolute_calls_scratch,
              "SEG_MAF" => {
                "@sample_name.segtab.txt" => :absolute_segs_scratch
              }
            }
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
          "@sample_name" => :qc_coverage_base,
          "fastqc" => {
            "@sample_name.@fastq_name.pdf" => :fastqc_pdf
          }
        },
        "@cohort_name.duplication_metrics" => :duplication_metrics,
        "@cohort_name.qc_summary" => :qc_summary
      },
      ":output_dir" => {
        "@sample_name" => {
          "@sample_name.:bam_label.bam" => :output_bam,
          "@sample_name.:bam_label.bam.bai" => :sample_bai,
          "@sample_name.mut.txt" => :tumor_muts,
          "@sample_name.stratton.pdf" => :stratton_plot_pdf,
          "@sample_name.somatic.maf" => :tumor_maf,
          "@sample_name.germline.maf" => :germline_maf,
          "@sample_name.all_muts.maf" => :all_muts_maf,
          "@sample_name.gene_cnr" => :tumor_gene_cnr,
          "@sample_name.exon_cnr" => :sample_exon_cnr,
          "@sample_name.cnvkit.cnr" => :sample_cnr,
          "@sample_name.cnr.RData" => :tumor_cnr_rdata,
          "@sample_name.cnvkit.RData" => :tumor_cnvkit_rdata,
          "@sample_name.cnvkit.seg" => :tumor_cnvkit_seg,
          "@sample_name.ascat.RData" => :tumor_ascat_rdata,
          "@sample_name.ascat.txt" => :tumor_ascat_txt,
          "@sample_name.cnr.seg" => :tumor_cnr_seg,
          "@sample_name.absolute.seg" => :tumor_absolute_seg,
          "@sample_name.ABSOLUTE.RData" => :absolute_rdata,
          "@sample_name.ABSOLUTE.plot.pdf" => :absolute_pdf,
          "@sample_name.mutations" => :tumor_mutations,
          "@sample_name.univ_geno_muts.vcf" => :ug_muts,
        },
        "@cohort_name" => {
          "@cohort_name.ABSOLUTE.table.txt" => :absolute_calls,
          "@cohort_name.seg" => :combined_seg,
          "@cohort_name.somatic.maf" => :combined_somatic_maf,
          "@cohort_name.absolute.pdf" => :combined_absolute_pdf,
          "@cohort_name.main_log" => :main_log_copy,
          "@cohort_name.qc_summary" => :qc_summary_copy
        },
        ":normal_name" => {
          ":normal_name.exon_cnr" => :normal_exon_cnr
        }
      }
    })

    def_var :tumor_somatic_mafs do tumor_samples.map{|s| tumor_maf(s) } end
    def_var :tumor_cnr_segs do tumor_samples.map{|s| tumor_cnr_seg(s) } end
    def_var :tumor_absolute_pdfs do tumor_samples.map{|s| absolute_pdf(s) } end

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
              i.fq1 ?  [ make_fastq(i.fq1), make_fastq(i.fq2) ] : nil
            }.flatten
        end
        s.extend_with :chroms => chromosomes
        s.extend_with :chunks => make_chunks(s)
        s.add_member :lane_name, "lane0" if !s.lane_name
        s.add_member :patient_name, "patient0" if !s.patient_name
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
      samples.group_by(&:lane_name).map do |lane,samples|
        { :lane_name => lane, :samples => samples }
      end
    end

    def make_patients
      samples.group_by(&:patient_name).map do |patient,samples|
        { :patient_name => patient, :samples => samples }
      end
    end

    def_var :lane_name do |l| (l||job_item).property :lane_name end
    def_var :chunk_size do 4_000_000 end
    def_var :disable_mutect_normal_filter do nil end
    job_items :chrom, :lane, :patient, :chunk, :fastq

    def_var :reads1_fastqs do |s| (s || sample).inputs.map{|input| input_fastq1(input) || reads1_fastq(input) } end
    def_var :reads2_fastqs do |s| (s || sample).inputs.map{|input| input_fastq2(input) || reads2_fastq(input) } end

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

    # Patient Merge
    def_var :patient_raw_sample_bams do patient.samples.map {|s| find_prop :patient_raw_sample_bam, s } end
    def_var :patient_raw_sample_bam do |s| find_prop :recal_bam, s end

    # PatientSplit
    def_var :realigned_patient_bams do patient.chroms.map{|c| realigned_patient_bam(c) } end
    def_var :patient_split_bams do patient.samples.map{|s| patient_split_bam s.chroms.find{|c| c.chrom_name == chrom.chrom_name} } end
    def_var :patient_sample_bams do patient.samples.map{|s| sample_bam(s) } end

    def_var :sample_patient_bams do sample.chroms.map{|c| patient_split_bam c } end

    # Hybrid qc
    def_var :qc_bam do tumor_bam end

    def_var :segment_smoothing do 1.5 end

    #mut_filter
    def_var :pindel_vcfs do sample.chroms.map{|c| pindel_vcf c } end
    def_var :mutect_snvses do sample.chroms.map{|c| mutect_snvs c } end

    def_var :normal_sample_covs do normal_samples.map{|s| sample_cov s } end

    # Absolute
    def_var :absolute_rdatas do tumor_samples.map{|s| absolute_rdata s } end

    def_var :mutations_config do "#{config_dir}/exome_mutations.yml" end
  end
end
