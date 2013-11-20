#!/bin/bash
module Exome
  class UnivGenoNormals
    include Pipeline::Step
    runs_tasks :unified_genotyper, :annotate_muts, :quality_filter, :filter_muts
    runs_on :samples
    resources :threads => 12

    class UnifiedGenotyper
      include Pipeline::Task
      requires_files :sample_bam, :interval_list
      outs_files :ug_raw_vcf

      def run
	log_info "Running Unified Genotyper"
	gatk :unified_genotyper,
		:genotype_likelihoods_model => :BOTH, 
                :genotyping_mode => :DISCOVERY,
		:input_file => config.sample_bam,
		:dbsnp => config.dbsnp_vcf,
		:intervals => config.interval_list, :baq => :CALCULATE_AS_NECESSARY,
		:standard_min_confidence_threshold_for_calling => 30.0,
		:standard_min_confidence_threshold_for_emitting => 10.0,
		:min_base_quality_score => 20, :output_mode => :EMIT_VARIANTS_ONLY,
		:out => config.ug_raw_vcf or error_exit "Unified Genotyper SNP calling failed"

      end
    end
    class AnnotateMuts
      include Pipeline::Task
      requires_files :sample_bam, :ug_raw_vcf
      dumps_file :ug_annotated_vcf

      def run
	log_info "Annotating Unified Genotyper SNPs"
	gatk :variant_annotator,
		:input_file => config.sample_bam,
                :num_threads => 1,
		:dbsnp => config.dbsnp_vcf,
		:intervals => config.ug_raw_vcf,
		:variant => config.ug_raw_vcf,
		:baq =>  :CALCULATE_AS_NECESSARY,
		:annotation => [ :QualByDepth, :RMSMappingQuality, :MappingQualityZero, :LowMQ, :MappingQualityRankSumTest, :FisherStrand, :HaplotypeScore, :ReadPosRankSumTest, :Coverage ],
		:out => config.ug_annotated_vcf or error_exit "Unified Genotyper SNP annotation failed"
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
                :num_threads => 1,
		:baq => :CALCULATE_AS_NECESSARY,
		:filterExpression => [ '"QD < 2.0"', '"MQ < 40.0"', '"FS > 60.0"', '"HaplotypeScore > 13.0"', '"MQRankSum < -12.5"', '"ReadPosRankSum < -8.0"' ],
		:filterName => [ :QDFilter, :MQFilter, :FSFilter, :HaplotypeScoreFilter, :MQRankSumFilter, :ReadPosFilter ],
		:out => config.ug_filtered_vcf or error_exit "Unified Genotyper SNP filtration failed"

      end
    end
    class FilterMuts
      include Pipeline::Task
      requires_files :ug_filtered_vcf
      outs_file :normal_muts

      def run
        muts = []
        headers = [ :gene, :sample, :chrom, :pos, :ref_allele, :alt_allele, :ref_count, :alt_count, :var_freq, :variant_classification, :protein_change, :transcript_change, :polyphen2_class, :cosmic_mutations, :dbSNP_RS ]
        v = VCF.read config.ug_filtered_vcf, config.mutations_config
        v.each do |l|
          next if l.skip_genotype?([:univ_geno_normal, :vcf] => config.sample_name)
          next if l.skip_oncotator?([:univ_geno_normal, :oncotator])
          muts.push :gene =>           l.onco.txp_gene,
            :sample => config.sample_name,
            :chrom =>                   l.chrom,
            :pos => l.pos,
            :ref_allele =>              l.ref,
            :alt_allele =>              l.alt,
            :var_freq => l.genotype(config.sample_name).alt_freq,
            :ref_count => l.genotype(config.sample_name).ref_count,
            :alt_count => l.genotype(config.sample_name).alt_count,
            :variant_classification =>  l.onco.txp_variant_classification,
            :protein_change =>          l.onco.txp_protein_change,
            :transcript_change => l.onco.txp_transcript_change, 
            :polyphen2_class  =>                  l.onco.pph2_class,
            :cosmic_mutations => l.onco.Cosmic_overlapping_mutations,
            :dbSNP_RS => l.onco.is_snp ? l.onco.dbSNP_RS : nil
        end
        File.open(config.normal_muts, "w") do |f|
          f.puts "#{headers.join("\t")}\n#{
            muts.sort_by{|m| -1 * m[:var_freq]}.map do |m|
              headers.map{|h| m[h] || "-" }.join("\t")
            end.join("\n")
          }"
        end
      end
    end
  end
end
