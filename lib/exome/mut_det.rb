#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

module Exome
  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :pindel, :pindel_vcf, :patch_pindel_vcf
    has_tasks :mutect, :pindel, :pindel_vcf, :patch_pindel_vcf, :strelka, :somatic_indel_detector, :patch_somatic_indel_vcf
    resources :threads => 1
    runs_on :tumor_samples, :chroms

    def vacuum
      config.sample_names.each do |s|
        FileUtils.rm Dir.glob("#{config.scratch_dir}/#{s}/AnnoVar*")
      end
    end

    class Mutect
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_files :mutect_snvs, :mutect_coverage

      def run
	      log_info "Running muTect for tumor #{config.sample_name}, normal #{config.normal_name}"
        mutect "input_file:normal" => config.normal_bam, "input_file:tumor" => config.tumor_bam,
          :normal_sample_name => config.normal_name,
          :tumor_sample_name => config.sample_name,
          :intervals => config.chrom.chrom_name,
          :no_normal_filter => config.disable_mutect_normal_filter,
          :out => config.mutect_snvs_tmp, :coverage_file => config.mutect_coverage or error_exit "muTect failed"

        # kludge to make sure mutect completes before ensuring this step
        FileUtils.mv config.mutect_snvs_tmp, config.mutect_snvs
      end
    end

    class Pindel
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_files :pindel_snv_d

      def run
        log_info "Running pindel"
        pindel :bams => [ 
            { :bam => config.tumor_bam, :name => config.sample_name},
            { :bam => config.normal_bam, :name => config.normal_name } ],
          :tempfile => config.pindel_list,
          :chromosome => config.chrom.chrom_name, :"output-prefix" => config.pindel_snvs or error_exit "Pindel failed"
      end
    end

    class PindelVcf
      include Pipeline::Task
      dumps_file :pindel_unpatched_vcf

      def run
        pindel_to_vcf :pindel_output_root => config.pindel_snvs, :vcf => config.pindel_unpatched_vcf or error_exit "Pindel2VCF failed"
      end
    end

    class PatchPindelVcf
      include Pipeline::Task
      requires_file :pindel_unpatched_vcf, :normal_bam, :tumor_bam
      dumps_file :pindel_vcf

      def run
        unpatched = VCF.new config.pindel_unpatched_vcf
        unpatched.preamble_lines.find{|i| i=~ /ID=AD/}.replace "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n"
        unpatched.each do |l|
          # skip invalid ones
          if l.alt == "<INS>"
            l.invalid = true
            next
          end
          l.format.push :DP 
          # get the actual depth for this locus
          normal_depth = count_depth config.normal_bam, l.seqname, l.pos
          tumor_depth = count_depth config.tumor_bam, l.seqname, l.pos

          support = l.genotype(config.normal_name).ad.to_i
          if support > normal_depth
            log_error "Mpileup depth is less than pindel supporting read count!"
            log_error "#{support} #{normal_depth} #{l.seqname}:#{l.pos}"
            normal_depth = support
          end
          l.genotype(config.normal_name).dp = normal_depth
          l.genotype(config.normal_name).ad = "#{normal_depth - support},#{support}"

          support = l.genotype(config.sample_name).ad.to_i
          if support > tumor_depth
            log_error "Mpileup depth is less than pindel supporting read count!"
            log_error "#{support} #{tumor_depth} #{l.seqname}:#{l.pos}"
            tumor_depth = support
          end
          l.genotype(config.sample_name).dp = tumor_depth
          l.genotype(config.sample_name).ad = "#{tumor_depth - support},#{support}"
        end
        unpatched.write config.pindel_vcf
      end
    end

    class SomaticIndelDetector
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      outs_file :somaticindel_unpatched_vcf
      outs_file :somaticindel_verbose => :exists

      def run
        somatic_indel_detector :"input_file:normal" => config.normal_bam,
          :"input_file:tumor" => config.tumor_bam,
          :verboseOutput => config.somaticindel_verbose,
          :out => config.somaticindel_unpatched_vcf,
          :filter_expressions => "'T_COV<6||N_COV<4||T_INDEL_CF<0.7'",
          :intervals => config.chrom.chrom_name or error_exit "SomaticIndelDetector failed"
      end
    end

    class PatchSomaticIndelVcf
      include Pipeline::Task
      requires_file :somaticindel_unpatched_vcf
      requires_file :somaticindel_verbose => :exists
      dumps_file :somaticindel_vcf

      def run
        verbose = SomaticIndelOut.new config.somaticindel_verbose
        unpatched = VCF.new
        unpatched.parse config.somaticindel_unpatched_vcf
        unpatched.each do |l|
          mut = verbose[ [ l.seqname, l.start, l.stop ] ]
          next if !mut
          l.genotype(config.normal_name.to_sym).dp = mut.normal_depth
          l.genotype(config.normal_name.to_sym).ad = mut.normal_allelic_depth

          l.genotype(config.sample_name.to_sym).dp = mut.tumor_depth
          l.genotype(config.sample_name.to_sym).ad = mut.tumor_allelic_depth
        end
        unpatched.write config.somaticindel_vcf
      end
    end

    class Strelka
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_file :strelka_vcf, :strelka_bam
      
      def run
        strelka :tumor => config.tumor_bam, :normal => config.normal_bam, :output => config.strelka_vcf, :chrom => config.chrom.chrom_name, :realigned_bam => config.strelka_bam or error_exit "Strelka failed"
      end
    end
  end
end
