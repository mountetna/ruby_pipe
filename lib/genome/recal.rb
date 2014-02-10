#!/usr/bin/env ruby
module Genome
  class LaneRecal
    include Pipeline::Step
    runs_tasks :count_covariates
    runs_on :lanes

    class CountCovariates
      include Pipeline::Task
      requires_file :lane_raw_sample_bams
      dumps_file :recal_grp

      def run
	log_info "Base-quality recalibration: Count covariates"
	gatk :base_recalibrator, :knownSites => config.reference_snp_vcf, 
		:input_file => config.lane_raw_sample_bams,
                :num_threads => nil,
		:out => config.recal_grp or error_exit "First CountCovariates failed"
      end
    end
  end

  class TableRecal
    include Pipeline::Step
    runs_tasks :table_recal
    runs_on :samples
    resources :threads => 1

    class TableRecal
      include Pipeline::Task
      requires_file :recal_grp, :raw_sample_bam
      dumps_file :recal_bam

      def run
	log_info "Base-quality recalibration: Table Recalibration"
	gatk :print_reads,
		:BQSR => config.recal_grp,
		:input_file => config.raw_sample_bam,
                :num_threads => nil,
		:out => config.recal_bam or error_exit "TableRecalibration failed"
      end
    end
  end

  class PatientRealign
    include Pipeline::Step
    runs_tasks :create_intervals, :realign_indels, :split_bam
    runs_on :patients, :chroms
    audit_report :patient_name
    resources :threads => 1

    class CreateIntervals
      include Pipeline::Task
      requires_file :patient_raw_sample_bams
      dumps_file :patient_intervals

      def run
	log_info "Creating intervals for indel detection"
	gatk :realigner_target_creator,
          :known => config.reference_indel_vcf,
          :intervals => config.chrom.chrom_name,
          :num_threads => nil,
          :input_file => config.patient_raw_sample_bams,
          :out => config.patient_intervals or error_exit "Interval creation failed"
      end
    end

    class RealignIndels
      include Pipeline::Task
      requires_file :patient_raw_sample_bams, :patient_intervals
      dumps_file :realigned_patient_bam

      def run
	log_info "Indel realignment"
	gatk :indel_realigner,
			:knownAlleles => config.reference_indel_vcf,
                        :consensusDeterminationModel => :USE_READS,
                        :intervals => config.chrom.chrom_name,
			:input_file => config.patient_raw_sample_bams,
                        :num_threads => nil,
			:targetIntervals => config.patient_intervals,
			:out => config.realigned_patient_bam or error_exit "Indel realignment failed"
      end
    end

    class SplitBam
      include Pipeline::Task
      requires_file :realigned_patient_bam
      dumps_file :patient_split_bams

      def run
	log_info "Split recalibrated bam into per-sample output bams"
	gatk :split_sam_file, :input_file => config.realigned_patient_bam, :outputRoot => config.patient_split_bam_root or error_exit "Splitting bam files failed"
      end
    end
  end

  class MakeSamples
    include Pipeline::Step
    runs_tasks :merge_files
    runs_on :samples

    class MergeFiles
      include Pipeline::Task
      requires_file :sample_patient_bams
      outs_file :sample_bam

      def run
        picard :merge_sam_files, :CREATE_INDEX => :true, :OUTPUT => config.sample_bam, :INPUT => config.sample_patient_bams or error_exit "bam file merge failed."
      end
    end
  end
end
