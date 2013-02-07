#!/usr/bin/env ruby
module Rna
  class UnivGeno
    include Pipeline::Step
    runs_tasks :unified_genotyper, :annotate_muts, :filter_muts

    class UnifiedGenotyper
      include Pipeline::Task
      requires_files :output_bams, :interval_list
      outs_files :ug_raw_vcf

      def run
	log_info "Running Unified Genotyper"
	gatk :unified_genotyper,
		:genotype_likelihoods_model => :BOTH, :genotyping_mode => :DISCOVERY,
		:input_file => config.output_bams,
		:dbsnp => config.dbsnp_vcf, :logging_level => :WARN,
		:intervals => config.interval_list, :baq => :CALCULATE_AS_NECESSARY,
		:standard_min_confidence_threshold_for_calling => 30.0,
		:standard_min_confidence_threshold_for_emitting => 10.0,
		:min_base_quality_score => 20, :output_mode => :EMIT_VARIANTS_ONLY,
		:out => config.ug_raw_vcf or error_exit "Unified Genotyper SNP calling failed"

      end
    end
    class AnnotateMuts
      include Pipeline::Task
      requires_files :output_bams, :ug_raw_vcf
      dumps_file :ug_annotated_vcf

      def run
	log_info "Annotating Unified Genotyper SNPs"
	gatk :variant_annotator,
		:input_file => config.output_bams,
		:dbsnp => config.dbsnp_vcf,
		:intervals => config.ug_raw_vcf,
		:variant => config.ug_raw_vcf,
		:baq =>  :CALCULATE_AS_NECESSARY,
		:annotation => [ :QualByDepth, :RMSMappingQuality, :MappingQualityZero, :LowMQ, :MappingQualityRankSumTest, :FisherStrand, :HaplotypeScore, :ReadPosRankSumTest, :DepthOfCoverage ],
		:out => config.ug_annotated_vcf or error_exit "Unified Genotyper SNP annotation failed"
      end
    end

    class FilterMuts
      include Pipeline::Task
      requires_files :ug_annotated_vcf
      dumps_file :ug_filtered_vcf

      def run
	log_info "Filtering Unified Genotyper SNPs"
	gatk :variant_filtration,
		:variant => config.ug_annotated_vcf,
		:baq => :CALCULATE_AS_NECESSARY,
		:filterExpression => [ '"QD < 2.0"', '"MQ < 40.0"', '"FS > 60.0"', '"HaplotypeScore > 13.0"', '"MQRankSum < -12.5"', '"ReadPosRankSum < -8.0"' ],
		:filterName => [ :QDFilter, :MQFilter, :FSFilter, :HaplotypeScoreFilter, :MQRankSumFilter, :ReadPosFilter ],
		:out => config.ug_filtered_vcf or error_exit "Unified Genotyper SNP filtration failed"

      end
    end
  end
end
