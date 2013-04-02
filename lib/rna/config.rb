module Rna
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def init_hook
      # add various bells and whistles here
      samples.each do |s|
        s.replicates.each do |r|
          r.add_member :replicate_name, "r#{r.index}"
        end
      end
    end

    def_var :sample_names do samples.collect(&:sample_name) end
    def_var :sample_bam do |s| input_bam(s) || output_bam(s) end
    def_var :sample_bams do samples.map{ |s| sample_bam(s) } end
    def_var :bam_label do "aligned.merged.sorted" end
    def_var :input_bam do |s| s ? s.input_bam : nil end
    def_var :replicate_bams do |s| (s||sample).replicates.map{ |r| replicate_bam r} end
    def_var :normal_bams do replicate_bams normal end
    def_var :replicate_bam do |r| input_bam(r) || output_bam(r) end
    def_var :normal_name do sample.normal_name || sample_names.first end
    def_var :normal do samples.find{|s| s.sample_name == normal_name} end
    def_var :tumor_bam do sample_bam(sample_name) end

    def_var :replicate do job_item end
    def_var :replicates do samples.collect(&:replicates).flatten end

    dir_tree({
      ":scratch_dir" => {
        "@sample_name" => {
          "@replicate_name" => {
            "tophat" => {
              "accepted_hits.bam" => :accepted_bam,
              "unmapped.bam" => :unmapped_bam,
              "merged.bam" => :merged_bam,
              "merged.sorted.bam" => :sorted_bam,
              "sorted_header.txt" => :sorted_header
            },
            "cufflinks" => {
              "transcripts.gtf" => :transcripts_gtf,
            }
          }
        },
        "@cohort_name" => {
          "cuffcompare" => {
            "cuffcomp" => :cuffcompare_prefix,
            "cuffcomp.tracking" => :tracking_file
          },
          "assembly.list" => :assembly_list,
          "merged.gtf" => :assembly_gtf
        }
      },
      ":metrics_dir" => {
        "@sample_name" => {
          "@sample_name.@replicate_name.flagstat" => :qc_flag,
          "@sample_name.@replicate_name.rnaseq_metrics" => :qc_rnaseq,
          "@sample_name.@replicate_name.rnaseq.pdf" => :qc_pdf,
          "@sample_name.@replicate_name.alignment_metrics" => :qc_align_metrics
        }
      },
      ":output_dir" => {
        "@sample_name" => {
          "@sample_name.@replicate_name.transcripts.gtf" => :output_gtf,
          "@sample_name.@replicate_name.:bam_label.bam" => :output_bam,
          "@sample_name.mutations" => :sample_mutations
        },
        "@cohort_name" => {
          "@cohort_name.ug.raw.vcf" => :ug_raw_vcf,
          "@cohort_name.ug.annotated.vcf" => :ug_annotated_vcf,
          "@cohort_name.ug.filtered.vcf" => :ug_filtered_vcf 
        }
      }
    })

    # tophat_align
    def_var :input_fastq1s do replicate.inputs.map(&:fq1) end
    def_var :input_fastq2s do replicate.inputs.map(&:fq2) end

    # count_transcripts
    def_var :transcripts_gtfs do replicates.map{ |r| transcripts_gtf r } end

    # qc
    def_var :qc_bam do replicate_bam end
    
    def_var :mutations_config do "#{config_dir}/rna_mutations.yml" end

    def_var :find_novel_transcripts do nil end
  end
end
