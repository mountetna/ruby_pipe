#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'
require 'seg'
require 'filter'

class TumorMaf < Maf
  EXTRA_HEADERS = [ :t_ref_count, :t_alt_count, :t_var_freq,
    :n_ref_count, :n_alt_count,
    :protein_change, :transcript_change,
    :segment_logr
  ]
  columns *(Maf.default_opts[:required] + EXTRA_HEADERS)
end

class AbsoluteMaf < Maf
  ABSOLUTE_HEADERS = [
    :t_ref_count, :t_alt_count, :t_var_freq,
    :n_ref_count, :n_alt_count, :protein_change,
    :transcript_change, :segment_logr
  ]
  columns *(Maf.default_opts[:required] + ABSOLUTE_HEADERS)
  display_names :start_position => :Start_position, :end_position => :End_position
end

module Exome
  class MutFilter
    include Pipeline::Step
    runs_tasks :filter_muts_pindel
    has_tasks :mutect_to_vcf, :snp_eff_annotate_mutect_vcf, :snp_eff_annotate_somatic_indel,
      :filter_muts_pindel, :concat_chroms, :filter_muts_annovar, :filter_muts_somatic_indel

    runs_on :tumor_samples, :chroms
    resources :internet => 1

    class FilterMuts
      include Pipeline::Task

      def vcf_to_maf mut
        seg = @segs.find do |seg|
          seg.overlaps? mut
        end
        {
          :hugo_symbol => mut.best_effect.gene_name,
          :entrez_gene_id => mut.best_effect.gene_id,
          :ncbi_build => 37,
          :chromosome => mut.chrom,
          :start_position => mut.start,
          :end_position => mut.stop,
          :reference_allele => mut.ref,
          :tumor_seq_allele1 => mut.alt,
          :tumor_seq_allele2 => nil,
          :tumor_validation_allele1 => nil,
          :tumor_validation_allele2 => nil,
          :match_norm_seq_allele1 => nil,
          :match_norm_seq_allele2 => nil,
          :match_norm_validation_allele1 => nil,
          :match_norm_validation_allele2 => nil,
          :verification_status => nil,
          :validation_status => nil,
          :mutation_status => nil,
          :sequencing_phase => nil,
          :sequence_source => nil,
          :validation_method => nil,
          :score => nil,
          :center => "cbc.ucsf.edu",
          :strand => "+",
          :variant_classification => mut.best_effect.annotation,
          :variant_type => nil,
          :dbsnp_rs => mut.id,
          :dbsnp_val_status => nil,
          :tumor_sample_barcode => config.sample_name,
          :matched_norm_sample_barcode => config.normal_name,
          :bam_file => config.sample_bam,
          :protein_change => mut.best_effect.hgvs_p,
          :transcript_change => mut.best_effect.hgvs_c,
          :t_ref_count => mut.t_ref_count,
          :t_alt_count => mut.t_alt_count,
          :t_var_freq => mut.t_var_freq,
          :n_ref_count => mut.n_ref_count,
          :n_alt_count => mut.n_alt_count,
          :segment_logr => seg ? seg.seg_mean.round(5) : nil
        }
      end

      def create_mafs
        @somatic_maf = TumorMaf.new
        @germline_maf = TumorMaf.new
        @all_muts_maf = AbsoluteMaf.new
      end

      def write_mafs
        @somatic_maf.sort_by! {|l| -l.t_var_freq }
        @germline_maf.sort_by! {|l| -l.t_var_freq }
        @somatic_maf.write config.tumor_chrom_maf
        @germline_maf.write config.germline_chrom_maf
        @all_muts_maf.write config.all_muts_chrom_maf
      end

      def load_mutect_vcf chrom
        m = MutectSnpeffVCF.new normal_name: config.normal_name.to_sym, tumor_name: config.sample_name.to_sym
        somatic_filter = Filter.new(@filters[ :mutect ][ :somatic ])
        germline_filter = Filter.new(@filters[ :mutect ][ :germline ])
        snpeff_filter = Filter.new(@filters[ :snpeff ])

        m.parse(config.mutect_annotated_vcf(chrom))
        m.each do |l|
          is_somatic = somatic_filter.passes?(l)
          #is_germline = germline_filter.passes?(l)
          next unless is_somatic# || is_germline

          # after this is any covered mutation - we want to remember all of these somewhere
          mut = vcf_to_maf(l)

          if snpeff_filter.passes?(l.best_effect)
            if is_somatic
              @somatic_maf << mut
            #else
              #@germline_maf << mut
            end
          end
          @all_muts_maf << mut if is_somatic
        end
      end

      def load_indel_vcf chrom
        indels = SomaticIndelSnpeffVCF.new.parse indel_vcf(chrom)
        indel_normal_filter = Filter.new(@filters[ indel_caller ][:normal])
        indel_tumor_filter = Filter.new(@filters[ indel_caller ][:tumor])
        snpeff_filter = Filter.new(@filters[ :snpeff ])
        normal = config.normal_name.to_sym
        tumor = config.sample_name.to_sym

        indels.each do |l|
          next if l.alt.include?("N") || l.ref.include?("N")
          next unless indel_normal_filter.passes?(l.genotype(normal)) && indel_tumor_filter.passes?(l.genotype(tumor))
          if !l.best_effect
            log_info "No best effect for #{l.range}, effects #{l.effects.count}"
            next
          end
          next unless snpeff_filter.passes?(l.best_effect)
          mut = vcf_to_maf l
          @somatic_maf << mut
        end
      end

      def run
        create_mafs

        @segs = Seg.new
        @segs.parse config.tumor_cnr_seg

        @filters = YAML.load(File.read(config.mutations_config))
        load_mutect_vcf config.chrom
        load_indel_vcf config.chrom

        write_mafs
      end
    end

    class FilterMutsPindel < FilterMuts
      class_init
      requires_files :pindel_vcf, :mutect_snvs, :tumor_cnr_seg
      dumps_file :tumor_chrom_maf, :germline_chrom_maf, :all_muts_chrom_maf

      def indel_caller; :pindel; end
      def indel_vcf chrom; config.pindel_vcf chrom; end
    end

    class MutectToVcf
      include Pipeline::Task

      # get this guy in VCF format so you can annotate it
      requires_file :mutect_snvs
      dumps_file :mutect_vcf

      def mutect_to_vcf snv_file, vcf_file
        mt = MuTect.new
        mt.parse snv_file
        v = VCF.new
        normal = config.normal_name.to_sym
        tumor = config.sample_name.to_sym
        v.add_columns :format, normal, tumor
        v.samples.concat [ normal, tumor ]
        mt.each do |m|
          v << {
                 :chrom => m.contig, :pos => m.position, :ref => m.ref_allele,
                 :info => "JM=#{m.judgement};CV=#{m.covered};MQ0=#{m.map_q0_reads}",
                 :alt => m.alt_allele, :qual => ".", :filter => ".", :id => ".",
                 :format => [ "DP", "AD" ],
                 normal => [ m.n_ref_count + m.n_alt_count, [ m.n_ref_count, m.n_alt_count ].join(",") ],
                 tumor => [ m.t_ref_count + m.t_alt_count, [ m.t_ref_count, m.t_alt_count ].join(",") ]
               }
        end
        v.write vcf_file
      end

      def run
        mutect_to_vcf config.mutect_snvs, config.mutect_vcf
      end
    end

    class SnpEffAnnotateMutectVcf
      include Pipeline::Task

      requires_file :mutect_vcf
      dumps_file :mutect_annotated_vcf

      def run
        snpeff config.mutect_vcf, config.mutect_annotated_vcf
      end
    end

    class SnpEffAnnotateSomaticIndel
      include Pipeline::Task

      requires_file :somaticindel_vcf
      dumps_file :somaticindel_annotated_vcf

      def run
        snpeff config.somaticindel_vcf, config.somaticindel_annotated_vcf
      end
    end

    class FilterMutsSomaticIndel < FilterMuts
      class_init
      requires_files :somaticindel_annotated_vcf, :mutect_annotated_vcf, :tumor_cnr_seg
      dumps_file :tumor_chrom_maf, :germline_chrom_maf, :all_muts_chrom_maf

      def indel_caller; :somaticindel; end
      def indel_vcf chrom; config.somaticindel_vcf chrom; end
    end

    class ConcatChroms
      include Pipeline::Task
      requires_files :chroms__pindel_vcfs, :chroms__mutect_snvs
      dumps_file :mutect_all_snvs, :pindel_all_vcf

      def run
        mutect = nil
        config.chroms__mutect_snvs.each do |mf|
          m = MuTect.new mf
          if mutect
            mutect.lines.concat m.lines
          else
            mutect = m
          end
        end
        mutect.write config.mutect_all_snvs if mutect
        vcf = nil
        config.chroms__pindel_vcfs.each do |vf|
          v = VCF.new vf
          v.lines.select! do |l|
            (l.ref.size == 1 || l.alt.size == 1) && l.alt !~ /[^ATGC]/  && l.ref !~ /[^ATGC]/
          end
          if vcf
            vcf.lines.concat v.lines
          else
            vcf = v
          end
        end
        vcf.write config.pindel_all_vcf if vcf
      end
    end

    class FilterMutsAnnovar
      include Pipeline::Task
      requires_files :pindel_all_vcf, :mutect_all_snvs
      outs_file :tumor_muts

      def run
        log_info "Filtering mutect and indel output..."
        filter_muts config.mutect_all_snvs, config.pindel_all_vcf, config.tumor_muts or error_exit "Filtering failed"
      end
    end
  end

  class CombineMuts
    include Pipeline::Step
    runs_on :tumor_samples
    runs_tasks :concat_mafs

    class ConcatMafs
      include Pipeline::Task
      requires_files :chroms__tumor_chrom_mafs, :chroms__germline_chrom_mafs, :chroms__all_muts_chrom_mafs
      outs_files :tumor_maf, :germline_maf, :all_muts_maf

      def write_mafs
        @somatic_maf.sort_by! {|l| -l.t_var_freq.to_f }
        @germline_maf.sort_by! {|l| -l.t_var_freq.to_f }
        @somatic_maf.write config.tumor_maf
        @germline_maf.write config.germline_maf
        @all_muts_maf.write config.all_muts_maf
      end

      def run
        @somatic_maf = TumorMaf.new
        @germline_maf = TumorMaf.new
        @all_muts_maf = AbsoluteMaf.new

        config.sample.chroms.each do |chrom|
          m = TumorMaf.new.parse config.tumor_chrom_maf(chrom)
          @somatic_maf.concat m

          m = TumorMaf.new.parse config.germline_chrom_maf(chrom)
          @germline_maf.concat m

          m = AbsoluteMaf.new.parse config.all_muts_chrom_maf(chrom)
          @all_muts_maf.concat m
        end
        write_mafs
      end
    end
  end

end
