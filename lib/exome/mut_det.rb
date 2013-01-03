#!/usr/bin/env ruby
module Exome
  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :somatic_indels, :annotate_indels, :filter_muts
    class Mutect
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :interval_list
      outs_files :mutect_snvs, :mutect_coverage

      def run
	log_info "Running muTect..."
        mutect "input_file:normal" => config.normal_bam, "input_file:tumor" => config.tumor_bam,
          :out => config.mutect_snvs, :coverage_file => config.mutect_coverage or error_exit "muTect failed"
      end
    end
    class SomaticIndels
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :interval_list
      outs_file :mutect_indels_raw

      def run
	log_info "Running Somatic Indel Detector..."
	gatk :somatic_indel_detector, "input_file:normal" => config.normal_bam,
		"input_file:tumor" => config.tumor_bam,
		"intervals" => config.interval_list,
                "maxNumberOfReads" => 10000,
                "window_size" => 225,
                "filter_expressions" => '"N_COV<8||T_COV<14||T_INDEL_F<0.1||T_INDEL_CF<0.7"',
		"out" => config.mutect_indels_raw or error_exit "Indel detection failed"
      end
    end
    class AnnotateIndels
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :mutect_indels_raw, :interval_list
      outs_file :mutect_indels_anno
      
      def run
	log_info "Annotating raw indel calls..."
	gatk :variant_annotator, 
		:variant => config.mutect_indels_raw,
		:intervals => config.interval_list,
		:"input_file:normal" => config.normal_bam,
		:"input_file:tumor" => config.tumor_bam,
		:dbsnp => config.dbsnp_vcf,
		:group => "StandardAnnotation",
		:out => config.mutect_indels_anno or error_exit "Indel annotation failed"
      end
    end
    class FilterMuts
      include Pipeline::Task
      requires_files :mutect_snvs, :mutect_indels_anno
      dumps_file :mutect_indels_temp
      outs_file :mutect_mutations

      def reorder_vcf(file,tumor_name,normal_name,out_file)
        header = []
        data = []

        File.foreach(file).each do |line|
          line.chomp!
          next if line =~ /^##/
          line = line.split("\t")
          if line.first == "#CHROM"
            header = line
          else
            data.push line
          end
        end

        tumor = header.index(tumor_name)
        normal = header.index(normal_name)
        header[9] = tumor_name
        header[10] = normal_name
        File.open(out_file,"w") do |f|
          data.each do |line|
            t = line[tumor]
            n = line[normal]
            line[9] = t
            line[10] = n
            f.puts line.join("\t")
          end
        end
      end

      def run
        normal_name = sam_sample_name config.normal_bam
        tumor_name = sam_sample_name config.tumor_bam

	log_info "Reordering indel vcf"
	reorder_vcf config.mutect_indels_anno, tumor_name, normal_name, config.mutect_indels_temp or error_exit "Reordering failed"

	log_info "Filtering mutect and indel output..."
	filter_muts config.mutect_snvs, config.mutect_indels_temp, config.mutect_mutations or error_exit "Filtering failed"
      end
    end
  end
end
