#!/bin/bash
module Exome
  class UnivGenoCall
    include Pipeline::Step
    runs_tasks :unified_genotyper, :annotate_muts, :quality_filter
    runs_on :patients
    resources :threads => 12

    class UnifiedGenotyper
      include Pipeline::Task
      requires_files :patient_sample_bams, :interval_list
      outs_files :ug_raw_vcf

      def run
	log_info "Running Unified Genotyper"
	gatk :unified_genotyper,
		:genotype_likelihoods_model => :BOTH, 
                :genotyping_mode => :DISCOVERY,
		:input_file => config.patient_sample_bams,
		:dbsnp => config.reference_snp_vcf,
		:intervals => config.interval_list, :baq => :CALCULATE_AS_NECESSARY,
		:standard_min_confidence_threshold_for_calling => 30.0,
		:standard_min_confidence_threshold_for_emitting => 10.0,
		:min_base_quality_score => 20, :output_mode => :EMIT_VARIANTS_ONLY,
		:out => config.ug_raw_vcf or error_exit "Unified Genotyper SNP calling failed"

      end
    end
    class AnnotateMuts
      include Pipeline::Task
      requires_files :patient_sample_bams, :ug_raw_vcf
      dumps_file :ug_annotated_vcf

      def run
	log_info "Annotating Unified Genotyper SNPs"
	gatk :variant_annotator,
		:input_file => config.patient_sample_bams,
                :num_threads => 1,
		:dbsnp => config.reference_snp_vcf,
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
      # requires_files :snp_annotated_vcf
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
  end
  class UnivGenoAnnotate
    include Pipeline::Step
    runs_task :filter_muts
    runs_on :samples
    class FilterMuts
      include Pipeline::Task
      requires_files :ug_filtered_vcf
      outs_file :ug_muts

      def run
        muts = []
        v = VCF.read config.ug_filtered_vcf, config.mutations_config
        v.samples.replace [ config.sample_name ]
        v.print config.ug_muts
      end
    end
  end
end
