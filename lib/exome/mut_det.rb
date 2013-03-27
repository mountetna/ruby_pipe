#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
module Exome
  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :pindel, :pindel_vcf
    resources :threads => 12
    job_list do config.tumor_samples.map{|s| s.extend_with( :chroms => config.chroms ) }.flatten end

    def vacuum
      config.sample_names.each do |s|
        FileUtils.rm Dir.glob("#{config.scratch_dir}/#{s}/AnnoVar*")
      end
    end

    class Mutect
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :interval_list
      dumps_files :mutect_snvs, :mutect_coverage

      def run
	log_info "Running muTect for tumor #{config.sample_name}, normal #{config.normal_name}"
        mutect "input_file:normal" => config.normal_bam, "input_file:tumor" => config.tumor_bam,
          :intervals => config.chrom,
          :out => config.mutect_snvs, :coverage_file => config.mutect_coverage or error_exit "muTect failed"
      end
    end

    class Pindel
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam, :interval_list
      dumps_files :pindel_snv_d

      def run
        log_info "Running pindel"
        pindel :bams => [ 
            { :bam => config.tumor_bam, :name => config.sample_name},
            { :bam => config.normal_bam, :name => config.normal_name } ],
          :tempfile => config.pindel_list,
          :chromosome => config.chrom, :output_prefix => config.pindel_snvs or error_exit "Pindel failed"
      end
    end
    class PindelVcf
      include Pipeline::Task
      requires_files :pindel_snv_d
      dumps_file :pindel_vcf

      def run
        pindel_to_vcf :pindel_output_root => config.pindel_snvs, :vcf => config.pindel_vcf or error_exit "Pindel2VCF failed"
      end
    end
  end
  class MutFilter
    include Pipeline::Step
    runs_tasks :filter_muts
    resources :threads => 12
    job_list do config.tumor_samples end

    class FilterMuts
      include Pipeline::Task
      requires_files :pindel_vcfs, :mutect_snvses
      outs_file :tumor_maf

      def run
        File.open(config.tumor_maf, "w") do |f|
          f.puts [ :gene, :chrom, :pos, :ref_allele, :alt_allele, :tumor_ref_count, :tumor_alt_count, :normal_ref_count, :normal_alt_count, :variant_classification, :protein_change, :class ].join("\t")
          config.chroms.each do |chrom|
            m = MuTect.new config.mutect_snvs(chrom[:contig]), config.mutations_config
            m.each do |l|
              next if l.skip_mutect?
              next if l.skip_oncotator?
              f.puts [ l.onco.txp_gene, l.chrom, l.pos, l.ref, l.alt, 
                l.tumor_ref_count, l.tumor_alt_count, l.normal_ref_count, l.normal_alt_count, 
                l.onco.txp_variant_classification, l.onco.txp_protein_change, 
                l.onco.pph2_class ].join("\t")
            end
            v = VCF.new config.pindel_vcf(chrom[:contig]), config.mutations_config
            v.each do |l|
              next if l.skip_genotype?([:pindel, :normal] => config.normal_name) || l.skip_genotype?([:pindel, :tumor] => config.sample_name)
              next if l.skip_oncotator?
              f.puts [ l.onco.txp_gene, l.chrom, l.pos, l.ref, l.alt, 
                "-", l.genotype(config.normal_name).approx_depth, "-", l.genotype(config.sample_name).approx_depth, 
                l.onco.txp_variant_classification, l.onco.txp_protein_change, 
                l.onco.pph2_class ].join("\t")
            end
          end
        end
      end
    end
  end
end
