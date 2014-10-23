#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

module Exome
  class MutFilter
    include Pipeline::Step
    runs_tasks :filter_muts_pindel
    has_tasks :filter_muts_pindel, :concat_chroms, :filter_muts_annovar, :filter_muts_somatic_indel
    runs_on :tumor_samples, :chroms
    resources :internet => 1

    class FilterMuts
      include Pipeline::Task

      EXTRA_HEADERS = [
        :tumor_ref_count, :tumor_alt_count, :tumor_var_freq,
        :normal_ref_count, :normal_alt_count,
        :protein_change, :transcript_change,
        :polyphen2_class, :cosmic_mutations, :segment_logr
      ]
      ABSOLUTE_HEADERS = [
        :t_ref_count, :t_alt_count, :tumor_var_freq,
        :normal_ref_count, :normal_alt_count, :protein_change,
        :transcript_change, :polyphen2_class, :cosmic_mutations, :segment_logr
      ]

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
        @somatic_maf = Maf.new
        @germline_maf = Maf.new
        @all_muts_maf = Maf.new

        @somatic_maf.header.concat EXTRA_HEADERS
        @germline_maf.header.concat EXTRA_HEADERS
        @all_muts_maf.header.concat ABSOLUTE_HEADERS
        @all_muts_maf.header.map! do |l|
          l.to_s =~ /_Position/ ? l.to_s.sub(/_Position/,"_position").to_sym : l
        end
      end

      def write_mafs
        @somatic_maf.sort_by! {|l| -l.tumor_var_freq }
        @germline_maf.sort_by! {|l| -l.tumor_var_freq }
        @somatic_maf.write config.tumor_chrom_maf
        @germline_maf.write config.germline_chrom_maf
        @all_muts_maf.write config.all_muts_chrom_maf
      end

      def load_mutect_snvs chrom
        MuTect.new(config.mutect_snvs(chrom), mutation_config: config.mutations_config).each do |l|
          next unless l.keep_somatic? || l.keep_germline?
          log_info "Annotating #{l.contig}:#{l.position}"
          mut = mutect_to_maf l
          unless l.skip_oncotator?
            @somatic_maf << mut if l.keep_somatic?
            @germline_maf << mut if l.keep_germline? && l.mut.onco.Cosmic_overlapping_mutations
          end
          @all_muts_maf << mut if l.keep_somatic?
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

        @segs = HashTable.new config.tumor_cnr_seg, :header => { :ID => :str, :Chromosome => :str, :Start => :int, :End => :int, :Num_Probes => :int, :Segment_Mean => :float }

        load_mutect_snvs config.chrom
        load_indel_snvs config.chrom

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

    class FilterMutsSomaticIndel < FilterMuts
      class_init
      requires_files :somaticindel_vcf, :mutect_snvs, :tumor_cnr_seg
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
        @somatic_maf.sort_by! {|l| -l.tumor_var_freq.to_f }
        @germline_maf.sort_by! {|l| -l.tumor_var_freq.to_f }
        @somatic_maf.write config.tumor_maf
        @germline_maf.write config.germline_maf
        @all_muts_maf.write config.all_muts_maf
      end

      def run
        @somatic_maf = Maf.new
        @germline_maf = Maf.new
        @all_muts_maf = Maf.new

        config.sample.chroms.each do |chrom|
          m = Maf.new config.tumor_chrom_maf(chrom)
          @somatic_maf.header = m.header
          @somatic_maf.concat m

          m = Maf.new config.germline_chrom_maf(chrom)
          @germline_maf.header = m.header
          @germline_maf.concat m

          m = Maf.new config.all_muts_chrom_maf(chrom)
          @all_muts_maf.header = m.header
          @all_muts_maf.concat m
        end
        write_mafs
      end
    end
  end

end
