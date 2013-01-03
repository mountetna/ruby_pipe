#!/bin/bash

# This script generates Universal Genotyper annotations from a NORMAL/BAM pair

#PBS -o $OUTFILE.log

. $LIB_DIR/base_opts $CONFIG

require_var MUTECT_MUTATIONS NORMAL_BAM TUMOR_BAM

log_info "Annotating mutations"

task_generate_bed() {
	log_info "Generate bed file from mutations..."
	l_annotate_bed_ug $MUTECT_MUTATIONS > $annovar_temp_bed || error_exit "Bed file creation failed"
}
run_task generate_bed "$annovar_temp_bed" "$MUTECT_MUTATIONS"

task_generate_UG() {
	log_info "Generate Unified Genotyper data..."
	p_gatk UnifiedGenotyper \
		--genotype_likelihoods_model SNP \
		--genotyping_mode DISCOVERY \
		--input_file $NORMAL_BAM \
		--input_file $TUMOR_BAM \
		--dbsnp $DBSNP_VCF \
		--logging_level WARN \
		--intervals $annovar_temp_bed \
		-baq CALCULATE_AS_NECESSARY \
		--noSLOD \
		--standard_min_confidence_threshold_for_calling 30.0 \
		--standard_min_confidence_threshold_for_emitting 10.0 \
		--min_base_quality_score 20 \
		--output_mode EMIT_VARIANTS_ONLY \
		--out $ug_raw_snps_vcf || error_exit "Unified Genotyper SNP calling failed"
}
run_task generate_UG "$ug_raw_snps_vcf" "$NORMAL_BAM $TUMOR_BAM $annovar_temp_bed"

task_annotate_UG() {
	log_info "Annotating Unified Genotyper SNPs..."
	p_gatk VariantAnnotator \
		--input_file $NORMAL_BAM \
		--input_file $TUMOR_BAM \
		--dbsnp $DBSNP_VCF \
		--intervals $ug_raw_snps_vcf \
		--variant $ug_raw_snps_vcf \
		-baq CALCULATE_AS_NECESSARY \
		--annotation QualByDepth \
		--annotation RMSMappingQuality \
		--annotation MappingQualityZero \
		--annotation LowMQ \
		--annotation MappingQualityRankSumTest \
		--annotation FisherStrand \
		--annotation HaplotypeScore \
		--annotation ReadPosRankSumTest \
		--annotation DepthOfCoverage \
		--out $ug_snps_annotated_vcf || error_exit "Unified Genotyper SNP annotation failed"
}
run_task annotate_UG "$ug_snps_annotated_vcf" "$ug_raw_snps_vcf $NORMAL_BAN $TUMOR_BAM"

task_filter_UG() {
	log_info "Filtering Unified Genotyper SNPs..."
	p_gatk VariantFiltration \
		--variant $ug_snps_annotated_vcf \
		-baq CALCULATE_AS_NECESSARY \
		--filterExpression "QD < 2.0" \
		--filterName QDFilter \
		--filterExpression "MQ < 40.0" \
		--filterName MQFilter \
		--filterExpression "FS > 60.0" \
		--filterName FSFilter \
		--filterExpression "HaplotypeScore > 13.0" \
		--filterName HaplotypeScoreFilter \
		--filterExpression "MQRankSum < -12.5" \
		--filterName MQRankSumFilter \
		--filterExpression "ReadPosRankSum < -8.0" \
		--filterName ReadPosFilter	\
		--out $ug_snps_filtered_vcf || error_exit "Unified Genotyper SNP filtration failed"
}
run_task filter_UG "$ug_snps_filtered_vcf" "$ug_snps_annotated_vcf"

task_add_ug_data() {
	log_info "Add Unified Genotyper data..."
	annotate_UG $MUTECT_MUTATIONS $ug_snps_filtered_vcf > $temp1_mutations || error_exit "Unified Genotyper annotation failed"
}
run_task add_ug_data "$temp1_mutations" "$MUTECT_MUTATIONS $ug_snps_filtered_vcf"

task_add_cosmic_data() {
	log_info "Add COSMIC data..."

	annotate_COSMIC $temp1_mutations $COSMIC_FILE > $temp2_mutations || error_exit "COSMIC annotation failed"

	rm -f ${patientID}.${prefix}.temp1.mutations

	echo "[Annotate] Identify kinase genes..."
	$PYTHON /home/jocostello/shared/LG3_Pipeline/scripts/annotation_KINASE.py \
		${patientID}.${prefix}.temp2.mutations \
		$KINASEDATA \
		> ${patientID}.${prefix}.temp3.mutations || { echo "Kinase gene annotation failed"; exit 1; }

	rm -f ${patientID}.${prefix}.temp2.mutations

	echo "[Annotate] Identify cancer genes..."
	$PYTHON /home/jocostello/shared/LG3_Pipeline/scripts/annotation_CANCER.py \
		${patientID}.${prefix}.temp3.mutations \
		$CANCERDATA \
		$CONVERT \
		> ${patientID}.${prefix}.annotated.mutations || { echo "Cancer gene annotation failed"; exit 1; }

	rm -f ${patientID}.${prefix}.temp3.mutations

	echo "[Annotate] Finished!"
