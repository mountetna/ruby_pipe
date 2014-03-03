#!/usr/bin/env ruby

require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

module Genome
  class VariantDet
    include Pipeline::Step
    runs_tasks :variant_caller, :variant_annot, :quality_filter
    resources :threads => 1, :walltime => 50
    runs_on :patients, :chroms

    class VariantCaller
      include Pipeline::Task
      requires_files :samples__sample_bams
      dumps_files :ug_raw_vcf
      def run
        log_info "Running GATK UnifiedGenotypeCaller on normal/tumor sample"
        gatk :unified_genotyper, 
             :genotype_likelihoods_model => :BOTH,
             :genotyping_mode => :DISCOVERY,
             :intervals => config.chrom.chrom_name,
             :input_file => config.samples__sample_bams, 
             :"dbsnp" => config.reference_snp_vcf,
             :"standard_min_confidence_threshold_for_calling" => 30.0,
             :"standard_min_confidence_threshold_for_emitting" => 10.0,
             :"min_base_quality_score" => 20,
             :"output_mode" => "EMIT_VARIANTS_ONLY",
             :"out" => config.ug_raw_vcf or error_exit "GATK UnifiedGenotypeCaller failed"
      end
    end

    class VariantAnnot
      include Pipeline::Task
      requires_file :ug_raw_vcf
      outs_files :ug_annotated_vcf 

      def run
        log_info "Running GATK VariantRecal on VCFs"
        gatk :variant_annotator, :input_file => config.normal_bam, 
             :intervals => config.chrom.chrom_name, 
             :variant => config.ug_raw_vcf, :num_threads => 1,
             :dbsnp => config.reference_snp_vcf, 
             :baq =>  :CALCULATE_AS_NECESSARY,
             :annotation => [ :QualByDepth, :RMSMappingQuality,
                              :MappingQualityZero, :LowMQ,
                              :MappingQualityRankSumTest, :FisherStrand,
                              :HaplotypeScore, :ReadPosRankSumTest, 
                              :Coverage ],
             :out => config.ug_annotated_vcf or error_exit "GATK VariantRecalibrator failed" 
      end
    end

    class QualityFilter
      include Pipeline::Task
      requires_files :ug_annotated_vcf
      dumps_file :ug_filtered_vcf

      def run
        log_info "Filtering Unified Genotyper SNPs"
        gatk :variant_filtration,
                :variant => config.ug_annotated_vcf,
                :intervals => config.chrom.chrom_name,
                :num_threads => 1,
                :baq => :CALCULATE_AS_NECESSARY,
                :filterExpression => [ '"QD < 2.0"',
                                       '"MQ < 40.0"',
                                       '"FS > 60.0"',
                                       '"HaplotypeScore > 13.0"',
                                       '"MQRankSum < -12.5"',
                                       '"ReadPosRankSum < -8.0"' ],
                :filterName => [ :QDFilter, :MQFilter, :FSFilter,
                                 :HaplotypeScoreFilter, :MQRankSumFilter,
                                 :ReadPosFilter ],
                :out => config.ug_filtered_vcf or error_exit "Unified Genotyper SNP filtration failed"

      end
    end
  end

  class MergeVariants
    include Pipeline::Step
    runs_tasks :merge_variants
    resources :threads => 12, :walltime => 50
    runs_on :patients

    class MergeVariants
      include Pipeline::Task
      requires_files :chroms__ug_filtered_vcfs
      dumps_file :ug_vcf

      def run
        log_info "Combining annotated, filatered VCFs"
        gatk :combine_variants, :variant => config.chroms__ug_filtered_vcfs, :num_threads => 1,
             :out => config.ug_vcf, :"genotypemergeoption" => "UNSORTED" or error_exit "Merging UnifiedGenotyper failed"
      end
    end
  end 

  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :pindel, :pindel_vcf, :patch_pindel_vcf
    resources :threads => 1, :walltime => 50
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
          #intervals option removed because running on whole genome
	  :intervals => config.chrom, 
          :out => config.mutect_snvs_tmp, :coverage_file => config.mutect_coverage or error_exit "muTect failed"

        # kludge to make sure mutect completes before ensuring this step
        FileUtils.mv config.mutect_snvs_tmp, config.mutect_snvs
      end
    end
    #CHANGMT: this is for BWA-aligned sequences
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
          :chromosome => config.chrom, :"output-prefix" => config.pindel_snvs or error_exit "Pindel failed"
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
        unpatched = VCF.read config.pindel_unpatched_vcf
        unpatched.preamble_lines.find{|i| i=~ /ID=AD/}.replace "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n"
        log_info "Trying to find genotyes for #{config.normal_name}, #{config.sample_name}"
        unpatched.each do |l|
          #skip invalid ones
          if l.alt == "<INS>"
            l.invalid = true 
            next
          end
          l.format.push :DP 
          # get the actual depth for this locus
          normal_depth = count_depth config.normal_bam, l.mutation[:chrom], l.mutation[:pos]
          tumor_depth = count_depth config.tumor_bam, l.mutation[:chrom], l.mutation[:pos]
          support = l.genotype(config.normal_name).info[:AD].to_i
          if support > normal_depth
            log_error "Mpileup depth is less than pindel supporting read count!"
            log_error "#{support} #{normal_depth} #{l.chrom}:#{l.pos}" 
            normal_depth = support
          end
          l.genotype(config.normal_name).info[:DP] = normal_depth
          l.genotype(config.normal_name).info[:AD] = "#{normal_depth - support},#{support}"

          support = l.genotype(config.sample_name).info[:AD].to_i
          if support > tumor_depth
            log_error "Mpileup depth is less than pindel supporting read count!"
            log_error "#{support} #{tumor_depth} #{l.chrom}:#{l.pos}"
            tumor_depth = support
          end
          l.genotype(config.sample_name).info[:DP] = tumor_depth
          l.genotype(config.sample_name).info[:AD] = "#{tumor_depth - support},#{support}"
        end
        unpatched.write config.pindel_vcf
      end
    end
  end

  class MutFilter
    include Pipeline::Step
    runs_tasks :filter_muts
    has_tasks :filter_muts, :concat_chroms, :filter_muts_annovar
    runs_on :tumor_samples

    class FilterMuts
      include Pipeline::Task
      requires_files :pindel_vcfs, :mutect_snvses, :tumor_cnr_seg
      outs_file :tumor_maf, :germline_maf
      dumps_file :all_muts_maf

      def run
        somatic_maf = Maf.new
        germline_maf = Maf.new
        all_muts_maf = Maf.new
        addl = [ :tumor_ref_count, :tumor_alt_count, :tumor_var_freq, :normal_ref_count, :normal_alt_count, :protein_change, :transcript_change, :polyphen2_class, :cosmic_mutations, :segment_logr ]
        somatic_maf.headers.concat addl
        germline_maf.headers.concat addl

        # stupid fixes for absolute
        all_muts_maf.headers.concat [ :t_ref_count, :t_alt_count, :tumor_var_freq, :normal_ref_count, :normal_alt_count, :protein_change, :transcript_change, :polyphen2_class, :cosmic_mutations, :segment_logr ]
        all_muts_maf.headers.map! { |l| l.to_s =~ /_Position/ ? l.to_s.sub(/_Position/,"_position").to_sym : l }

        segs = HashTable.new config.tumor_cnr_seg
        config.sample.chroms.each do |chrom|
          MuTect.read(config.mutect_snvs(chrom), config.mutations_config).each do |l|
            next unless l.keep_somatic? || l.keep_germline?
            log_info "Annotating #{l.contig}:#{l.position}"
            seg = segs.find{|seg| seg[:Chromosome] == l.contig && seg[:Start].to_i < l.position.to_i && seg[:End].to_i > l.position.to_i}
            mut = {
              :hugo_symbol => l.onco.txp_gene, :center => "taylorlab.ucsf.edu",
              :ncbi_build => 37, :chromosome => l.contig.sub(/^chr/,""),
              :start_position => l.position, :end_position => l.position, :strand => "+",
              :variant_classification => l.onco.txp_variant_classification, :variant_type => l.onco.variant_type,
              :reference_allele => l.ref_allele, :tumor_seq_allele1 => l.alt_allele,
              :dbsnp_rs => (l.onco.is_snp ? l.onco.dbSNP_RS : nil), :dbsnp_val_status => l.onco.dbSNP_Val_Status,
              :tumor_sample_barcode => config.sample_name, :matched_norm_sample_barcode => config.normal_name,
              :bam_file => config.sample_bam,
              :tumor_ref_count => l.t_ref_count,
              :tumor_alt_count => l.t_alt_count,
              :t_ref_count => l.t_ref_count,
              :t_alt_count => l.t_alt_count,
              :tumor_var_freq => l.t_var_freq,
              :normal_ref_count => l.n_ref_count,
              :normal_alt_count => l.n_alt_count,
              :protein_change => l.onco.txp_protein_change,
              :transcript_change => l.onco.txp_transcript_change,
              :polyphen2_class => l.onco.pph2_class,
              :cosmic_mutations => l.onco.Cosmic_overlapping_mutations,
              :segment_logr => seg ? seg[:Segment_Mean].to_f.round(5) : nil
            }
            unless l.skip_oncotator?
              somatic_maf.add_line(mut) if l.keep_somatic?
              germline_maf.add_line(mut) if l.keep_germline? && l.onco.Cosmic_overlapping_mutations
            end
            all_muts_maf.add_line mut if l.keep_somatic?
          end

          v = VCF.read config.pindel_vcf(chrom), config.mutations_config
          v.each do |l|
	    log_info "Checking #{l.chrom}:#{l.pos}-#{l.end_pos}"
	    next if l.alt.include? "N"
            next if l.skip_genotype?([:pindel, :normal] => config.normal_name) || l.skip_genotype?([:pindel, :tumor] => config.sample_name)
            next if l.skip_oncotator?
            log_info "Annotating #{l.chrom}:#{l.pos}-#{l.end_pos}"
            seg = segs.find{|seg| seg[:Chromosome] == l.chrom && seg[:Start].to_i < l.pos.to_i && seg[:End].to_i > l.pos.to_i}
            mut = {
              :hugo_symbol => l.onco.txp_gene, :center => "taylorlab.ucsf.edu",
              :ncbi_build => 37, :chromosome => l.chrom.sub(/^chr/,""),
              :start_position => l.pos, :end_position => l.end_pos, :strand => "+",
              :variant_classification => l.onco.txp_variant_classification, :variant_type => l.onco.variant_type,
              :reference_allele => l.ref, :tumor_seq_allele1 => l.alt,
              :dbsnp_rs => (l.onco.is_snp ? l.onco.dbSNP_RS : nil), :dbsnp_val_status => l.onco.dbSNP_Val_Status,
              :tumor_sample_barcode => config.sample_name, :matched_norm_sample_barcode => config.normal_name,
              :bam_file => config.sample_bam,
              :tumor_ref_count => l.genotype(config.sample_name).ref_count,
              :tumor_alt_count => l.genotype(config.sample_name).alt_count,
              :t_ref_count => l.genotype(config.sample_name).ref_count,
              :t_alt_count => l.genotype(config.sample_name).alt_count,
              :tumor_var_freq => l.genotype(config.sample_name).alt_freq,
              :normal_ref_count => l.genotype(config.normal_name).ref_count,
              :normal_alt_count => l.genotype(config.normal_name).alt_count,
              :protein_change =>          l.onco.txp_protein_change,
              :transcript_change => l.onco.txp_transcript_change,
              :polyphen2_class  =>                  l.onco.pph2_class,
              :cosmic_mutations => l.onco.Cosmic_overlapping_mutations,
              :segment_logr => seg ? seg[:Segment_Mean].to_f.round(5) : nil
            }
            somatic_maf.add_line mut
          end
        end
        somatic_maf.sort_by! {|l| -l.tumor_var_freq }
        germline_maf.sort_by! {|l| -l.tumor_var_freq }
        somatic_maf.write config.tumor_maf
        germline_maf.write config.germline_maf
        all_muts_maf.write config.all_muts_maf
      end
    end

    class ConcatChroms
      include Pipeline::Task
      requires_files :pindel_vcfs, :mutect_snvses
      dumps_file :mutect_all_snvs, :pindel_all_vcf

      def run
        mutect = nil
        config.mutect_snvses.each do |mf|
          m = MuTect.read mf
          if mutect
            mutect.lines.concat m.lines
          else
            mutect = m
          end
        end
        mutect.write config.mutect_all_snvs if mutect
        vcf = nil
        config.pindel_vcfs.each do |vf|
          v = VCF.read vf
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
end
