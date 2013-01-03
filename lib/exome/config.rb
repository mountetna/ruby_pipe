module Exome
  class Config
    include Pipeline::Config

    def platform; "Illumina"; end
    def platform_unit; "Exome"; end

    # Align stuff
    def load_procs
      super

      @procs.update({
        :input_fastq1 => proc { input_fastq1_list.flatten[job_index] },
        :input_fastq1_list => proc { samples.map{ |s| s[:input_fastq1_list] }.flatten },
        :input_fastq2 => proc { input_fastq2_list.flatten[job_index] },
        :input_fastq2_list => proc { samples.map{ |s| s[:input_fastq2_list] }.flatten },
        :sample_tag => proc { samples.map{|s| s[:input_fastq1_list].to_enum.with_index.map{ |f,i| i } }.flatten[job_index] },
        :mated_sam => proc { "#{scratch}/#{sample_name}.#{sample_tag}.mateFixed.bwa.sam" },
        :read1_sai => proc { "#{scratch}/#{sample_name}.#{sample_tag}.read1.sai" },
        :read2_sai => proc { "#{scratch}/#{sample_name}.#{sample_tag}.read2.sai" },
        :paired_sam => proc { "#{scratch}/#{sample_name}.#{sample_tag}.paired.bwa.sam" },
        :renamed_sam => proc { "#{scratch}/#{sample_name}.#{sample_tag}.renamed.bwa.sam" },
        :renamed_bam => proc { "#{scratch}/#{sample_name}.#{sample_tag}.renamed.bwa.bam" },
        :sample_names => proc { samples.map{|s| s[:sample_name]}.flatten },
        :normal_name => proc { sample_names[0] },

        :sample_bam => proc { |s| "#{output_dir}/#{s || sample_name}/#{s || sample_name}.bwa.realigned.rmDups.recal.bam" },
        :sample_bams => proc { sample_names.map{|s| sample_bam(s) } },

        :output_bam => proc { |s| "#{output_dir}/#{s || sample_name}/#{s || sample_name}.bwa.realigned.rmDups.recal.bam" },
        :output_bams => proc { sample_names.map{|s| output_bam(s) } },

        # Merge stuff
        :aligned_bams => proc { samples[job_index][:input_fastq1_list].size.times.map { |t| "#{scratch}/#{sample_name}.#{t}.renamed.bwa.bam" } },
        :merged_bam => proc { "#{scratch}/merged.bam" },
        :merged_bai => proc { "#{scratch}/merged.bam.bai" },
          
        # Recal stuff
        :dedup_bam => proc { "#{scratch}/dedup.bam" },
        :mated_bam => proc { "#{scratch}/mated.bam" },
        :merged_intervals => proc { "#{scratch}/merged.intervals" },
        :realigned_bam => proc { "#{scratch}/realigned.bam" },
        :recal_bam => proc { "#{scratch}/recal.bam" },
        :recal_bai => proc { "#{scratch}/recal.bam.bai" },
        :recal_bams => proc { sample_bams },
        :recal_grp => proc { "#{scratch}/recal.grp" },
        :recal_metrics => proc { "#{metrics_dir}/#{sample_name}.recal_metrics" },
        :split_bam => proc { |s| s ? "#{scratch}/_splitbam_#{s}.bam" : "#{scratch}/_splitbam_" },
        :split_bams => proc { sample_names.map{ |s| split_bam(s) } },

        # hybrid_qc
        :qc_bam => proc { tumor_bam },
        :qc_flag => proc { "#{metrics_dir}/#{sample_name}.flagstat" },
        :qc_histogram => proc { "#{metrics_dir}/#{sample_name}.histogram" },
        :qc_hybrid => proc { "#{metrics_dir}/#{sample_name}.hybrid_selection_metrics" },
        :qc_inserts => proc { "#{metrics_dir}/#{sample_name}.insert_sizes" },

        # mut_det
        :mutect_snvs => proc { "#{output_dir}/#{sample_name}/#{sample_name}.snvs.raw.mutect.txt" },
        :mutect_indels_anno => proc { "#{output_dir}/#{sample_name}/#{sample_name}.indels.annotated.vcf" },
        :mutect_indels_raw => proc { "#{output_dir}/#{sample_name}/#{sample_name}.indels.raw.vcf" },
        :mutect_indels_temp => proc { "#{scratch}/#{sample_name}.indels.temp.vcf" },
        :mutect_mutations => proc { "#{output_dir}/#{sample_name}/#{sample_name}.mutations" },
        :mutect_coverage => proc { "#{output_dir}/#{sample_name}/#{sample_name}.snvs.coverage.mutect.wig" },

        :normal_bam => proc { sample_bam(normal_name) },
        :tumor_bam => proc { sample_bam(sample_name) },

        # copy_number
        :normal_cov => proc { },
        :tumor_cov => proc { }
      })
    end
  end
end
