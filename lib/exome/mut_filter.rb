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

      def mut_to_maf mut
        seg = @segs.find do |seg|
          seg.Chromosome.sub(/^chr/,'') == mut.short_chrom && seg.Start < mut.start && seg.End > mut.stop
        end
        {
          :hugo_symbol => mut.mut.onco.txp_gene,
          :entrez_gene_id => "",
          :ncbi_build => 37,
          :chromosome => mut.short_chrom,
          :start_position => mut.start,
          :end_position => mut.stop,
          :reference_allele => mut.mut.ref,
          :tumor_seq_allele1 => mut.mut.alt,
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
          :center => "taylorlab.ucsf.edu",
          :strand => "+",
          :variant_classification => mut.mut.onco.txp_variant_classification,
          :variant_type => mut.mut.onco.variant_type,
          :dbsnp_rs => (mut.mut.onco.is_snp ? mut.mut.onco.dbSNP_RS : nil),
          :dbsnp_val_status => mut.mut.onco.dbSNP_Val_Status,
          :tumor_sample_barcode => config.sample_name,
          :matched_norm_sample_barcode => config.normal_name,
          :bam_file => config.sample_bam,
          :protein_change => mut.mut.onco.txp_protein_change,
          :transcript_change => mut.mut.onco.txp_transcript_change,
          :polyphen2_class => mut.mut.onco.pph2_class,
          :cosmic_mutations => mut.mut.onco.Cosmic_overlapping_mutations,
          :segment_logr => seg ? seg.Segment_Mean.round(5) : nil
        }
      end

      def mutect_to_maf mut
        mut_to_maf(mut).merge({
          :tumor_ref_count => mut.t_ref_count,
          :tumor_alt_count => mut.t_alt_count,
          :t_ref_count => mut.t_ref_count,
          :t_alt_count => mut.t_alt_count,
          :tumor_var_freq => mut.t_var_freq,
          :normal_ref_count => mut.n_ref_count,
          :normal_alt_count => mut.n_alt_count,
        })
      end

      def indel_vcf_to_maf mut
        mut_to_maf(mut).merge({
          :tumor_ref_count => mut.genotype(config.sample_name).ref_count,
          :tumor_alt_count => mut.genotype(config.sample_name).alt_count,
          :t_ref_count => mut.genotype(config.sample_name).ref_count,
          :t_alt_count => mut.genotype(config.sample_name).alt_count,
          :tumor_var_freq => mut.genotype(config.sample_name).alt_freq,
          :normal_ref_count => mut.genotype(config.normal_name).ref_count,
          :normal_alt_count => mut.genotype(config.normal_name).alt_count,
        })
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
        mutation_filters = YAML.load(File.read(config.mutations_config))
        somatic_filter = Filter.new(mutation_filters[ :mutect ][ :somatic ])
        germline_filter = Filter.new(mutation_filters[ :mutect ][ :germline ])
        snpeff_filter = Filter.new(mutation_filters[ :snpeff ])

        m.parse(config.mutect_annotated_vcf(chrom))
        m.each do |l|
          is_somatic = somatic_filter.passes?(l)
          #is_germline = germline_filter.passes?(l)
          next unless is_somatic# || is_germline
          # after this is any covered mutation - we want to remember all of these somewhere
          
          mut = vcf_to_maf(l)

          if snpeff_filter.passes?(l.best_effect)
            # it has an interesting annotation
            if is_somatic
              @somatic_maf << mut
            #else
              #@germline_maf << mut
            end
          end
          # put all somatic mutations in this MAF for ABSOLUTE to use
          @all_muts_maf << mut if is_somatic
        end
      end

      def load_indel_snvs chrom
        VCF.new(indel_vcf(chrom), mutation_config: config.mutations_config).each do |l|
          begin
            log_info "Checking #{l.range}"
            next if l.alt.include?("N") || l.ref.include?("N")
            next if l.skip_genotype?([indel_caller, :normal] => config.normal_name) || l.skip_genotype?([indel_caller, :tumor] => config.sample_name)
            next if l.skip_oncotator?
            log_info "Annotating #{l.range}"
            mut = indel_vcf_to_maf l
            @somatic_maf << mut
          rescue ArgumentError => e
            log_info e.message
          end
        end
      end

      def run
        create_mafs

        @segs = Seg.new
        @segs.parse config.tumor_cnr_seg

        load_mutect_vcf config.chrom
        #load_indel_vcf config.chrom

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
                 :format => [ "AD", "DP" ],
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
