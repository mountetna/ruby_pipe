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
          "@{fastq_name}_fastqc" => {
            "." => :fastqc_output_dir,
            "summary.txt" => :fastqc_summary,
            "fastqc_report.html" => :fastqc_html
          },
          "input.@input_name.chaste.bam" => :chaste_bam,
          "input.@input_name.reads1.fastq.gz" => :reads1_fastq,
          "input.@input_name.reads2.fastq.gz" => :reads2_fastq,

          "input.@input_name.chunk_info" => :chunk_info,

          "input.read1.sai" => :read1_sai,
          "input.read2.sai" => :read2_sai,
          "@input_name.paired.bwa.sam" => :paired_sam,

          "chunk.@input_name.@chunk_name.read1.fastq.gz" => :chunk1_fastq,
          "chunk.@input_name.@chunk_name.read2.fastq.gz" => :chunk2_fastq,
          "chunk.@input_name.@chunk_name.read1.sai" => :read1_sai,
          "chunk.@input_name.@chunk_name.read2.sai" => :read2_sai,
          "chunk.@input_name.@chunk_name.mate_fixed.bwa.bam" => :mated_bam,
          "chunk.@input_name.@chunk_name.paired.bwa.sam" => :paired_sam,
          "chunk.@input_name.@chunk_name.dedup.bwa.bam" => :dedup_bam,
          "chunk.@input_name.@chunk_name.aligned.bwa.bam" => :aligned_bam,

          "@sample_name.@input_name.recal.bam" => :recal_bam,

          "@chrom_name.normal.snp.txt" => :normal_mut,

          "@sample_name.snps.baf"=> :tumor_baf,
          ":normal_name.snps.baf"=> :normal_baf,
          "@chrom_name.snvs.raw.mutect.txt" => :mutect_snv,
          "@chrom_name.snvs.raw.mutect.tmp.txt" => :mutect_snv_tmp,
          "@sample_name.snvs.raw.mutect.txt" => :mutect_all_snvs,
          "@chrom_name.snvs.coverage.mutect.wig" => :mutect_coverage,
          "@chrom_name.insert_mutations" => :insert_mutations,

          "@chrom_name.snvs.pindel" => :pindel_snvs,
          "@chrom_name.snvs.pindel_D" => :pindel_snv_d,
          "@chrom_name.pindel.conf" => :pindel_list,
          "@chrom_name.indels.unpatched.pindel.vcf" => :pindel_unpatched_vcf,
          "@chrom_name.indels.raw.pindel.vcf" => :pindel_vcf,
          "@chrom_name.indelocator.vcf" => :indelocator_unpatched_vcf,
          "@chrom_name.indels.raw.indelocator.vcf" => :indelocator_vcf,
          "@chrom_name.indelocator.txt" => :indelocator_output,

          "@chrom_name.somaticindel.unpatched.vcf" => :somaticindel_unpatched_vcf,
          "@chrom_name.indels.raw.somaticindel.vcf" => :somaticindel_vcf,
          "@chrom_name.somaticindel.verbose.txt" => :somaticindel_verbose,
          "@sample_name.indels.raw2.pindel.vcf" => :pindel_all_vcf,
          "@sample_name.ratio" => :tumor_ratio,
          "@chrom_name.somatic.maf" => :tumor_chrom_maf,
          "@chrom_name.germline.maf" => :germline_chrom_maf,
          "@chrom_name.all_muts.maf" => :all_muts_chrom_maf,
          "@sample_name.all_muts.maf" => :all_muts_maf,
          "@sample_name.pindel.output.txt" => :output_for_pindel,
          "@sample_name.temp.pindel.output.txt" => :temp_output_for_pindel,

          "@sample_name.@chrom_name.sclip.txt" => :sclip_file,
          "@sample_name.@chrom_name.cover" => :cover_file,
          "@sample_name.cover" => :sample_cover_file,
          "@sample_name.@chrom_name.predSV.txt" => :crest_rearr,
          "@sample_name.predSV.txt" => :sample_crest_rearr,

          ":normal_name.cov" => :normal_cov,
          "@sample_name.cov" => :tumor_cov,
          ":normal_name.gc_corrected.cov" => :normal_cov_gc,
          "@sample_name.gc_corrected.cov" => :tumor_cov_gc,
          "absolute" => {
            "." => :absolute_scratch,
            "@sample_name.ABSOLUTE.RData" => :absolute_rdata
          }
        },
        "@lane_name" => {
          "recal.grp" => :recal_grp,
          "recal_plot.pdf" => :recal_plot_pdf,
          "_splitbam_" => :lane_split_bam_root,
          "_splitbam_SAMPLE.bam" => :lane_split_bam
        },
        "@patient_name" => {
          "duplication_metrics" => :duplication_metrics,
          "@chrom_name.patient.intervals" => :patient_intervals,
          "@chrom_name.realigned_patient.bam" => :realigned_patient_bam,
          "@chrom_name.split." => :patient_split_bam_root,
          "@chrom_name.split.@sample_name.bam" => :patient_split_bam,
          "@patient_name.@chrom_name.ug.raw.vcf" => :ug_raw_vcf,
          "@patient_name.@chrom_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@patient_name.@chrom_name.ug.filtered.vcf" => :ug_filtered_vcf
        },
        "@cohort_name" => {
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
          "fastqc" => {
            "@sample_name.@fastq_name.pdf" => :fastqc_pdf
          }
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
          "@sample_name.gene_cnr" => :tumor_gene_cnr,
          "@sample_name.exon_cnr" => :tumor_exon_cnr,
          "@{sample_name}.cnr.RData" => :tumor_cnr_rdata,
          "@{sample_name}_pscbs.Rdata" => :tumor_cnr_rdata_pscbs,
          "@sample_name.pscbs.cnr.seg" => :tumor_cnr_seg_pscbs,
          "@sample_name.cnr.seg" => :tumor_cnr_seg,
          "@sample_name.mutations" => :tumor_mutations,
        },
        "@patient_name" => {
          "@patient_name.ug.filtered.vcf" => :ug_vcf
        }
      },
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
      samples.map(&:inputs).flatten.group_by(&:lane_name).map do |lane,inputs|
        { :lane_name => lane, :inputs => inputs }
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
            i.add_member :lane_name, "lane0" if !i.lane_name
            i.extend_with :chunks => make_chunks(i)
          end
          s.extend_with :fastqs => s.inputs.map{|i|
            [ make_fastq(i.fq1), make_fastq(i.fq2) ]
          }.flatten
        end
        s.extend_with :chroms => chromosomes
        s.add_member :patient_name, "patient#{s.patient || 0}"
      end
      @config.extend_with :lanes => make_lanes
      @config.extend_with :patients => make_patients
      patients.each do |patient|
        patient.extend_with :chroms => chromosomes
      end
      @config.extend_with :chroms => chromosomes
    end

    def_var :chunk_size do 4_000_000 end

    job_items :chrom, :lane, :patient, :chunk, :fastq, :input

    def_var :reads_bam do job_item.reads_bam end

    # Align
    def_var :input_fastq1 do |i| (i || job_item).property(:fq1) || (i || job_item).property(:reads1_fastq) end
    def_var :input_fastq2 do |i| (i || job_item).property(:fq2) || (i || job_item).property(:reads2_fastq) end

    # Lane Merge
    def_var :lane_aligned_bams do lane.inputs.map{|i| i.chunks.map{|c| aligned_bam c} }.flatten end

    # Patient Merge
    def_var :patient_recal_bams do
      patient.samples.map do |s|
        s.inputs.map do |i|
          recal_bam i
        end
      end.flatten
    end

    # PatientSplit
    def_var :patient_split_bams do patient.samples.map{|s| patient_split_bam s.chroms.find{|c| c.chrom_name == chrom.chrom_name} } end

    #mut_filter
    def_var :mutations_config do "#{config_dir}/genome_mutations.yml" end
  end
end
