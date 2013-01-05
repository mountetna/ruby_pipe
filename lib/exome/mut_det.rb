#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
module Exome
  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :somatic_indels, :annotate_indels, :filter_muts, :flag_inserts

    def vacuum
      config.sample_names.each do |s|
        FileUtils.rm Dir.glob("#{config.scratch_dir}/#{s}/AnnoVar*")
      end
    end

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
      dumps_file :mutect_indels_temp, :mutect_mutations

      def reorder_vcf(file,tumor_name,normal_name,out_file)
        data = HashTable.new(file,nil,"##")
        data.header[9] = tumor_name.to_sym
        data.header[10] = normal_name.to_sym
        data.print(out_file)
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
    class FlagInserts
      include Pipeline::Task
      requires_files :mutect_mutations, :qc_inserts, :tumor_bam
      dumps_file :insert_mutations

      def count_bad_inserts(chr, pos, low, high)
        reads = sam_reads(config.tumor_bam, chr, pos, pos)
        reads.count do |r|
          r[:insert_size].to_i.abs < low || r[:insert_size].to_i.abs > high || (r[:mate_contig] != '=' && r[:contig] != r[:mate_contig])
        end
      end

      def run
        # calculate insert size 5 and 95 percent intervals
        inserts = File.foreach(config.qc_inserts).select{ |s| s =~ /^[0-9]+\t[0-9]+$/ }.map{|s| s.chomp.split(/\t/).map(&:to_i)}

        total = inserts.inject(0) { |o,a| o += a.last }
        cumu_hist = inserts.each_with_object([0]).map { |a,o| [ a.first, o[0] += a.last] }
        i05 = cumu_hist.find{ |a| a.last > 0.05 * total }.first
        i95 = cumu_hist.find{ |a| a.last > 0.95 * total }.first

        # read the mutations into a damn data structure
        muts = HashTable.new(config.mutect_mutations)
        muts.header.insert( muts.header.index(:n_ref_count), :t_bad_inserts)

        muts.each do |l|
          l[:t_bad_inserts] = count_bad_inserts(l[:contig], l[:position], i05, i95)
        end

        muts.print(config.insert_mutations)
      end
    end
  end
end
