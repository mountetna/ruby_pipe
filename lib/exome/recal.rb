#!/usr/bin/env ruby
module Exome
  # the Broad folks do:
  # 1. Merge all BAMs for a lane (this might not be possible, but do your best) 
  # 2. recalibrate per lane
  # 3. Split and merge bams for each patient
  # 4. Dedup and indel realign
  # 5. Split and write bams
  class LaneMerge
    include Pipeline::Step
    audit_report :lane_name
    runs_tasks :merge_bam
    runs_on :lanes
    resources :threads => 1

    class MergeBam
      include Pipeline::Task
      requires_files :lane_raw_sample_bams
      dumps_file :merged_lane_bam

      def run
	log_info "Merging per-sample bam files for bulk recalibration"
        picard :merge_sam_files, :CREATE_INDEX => :true, :OUTPUT => config.merged_lane_bam, :INPUT => config.lane_raw_sample_bams or error_exit "bam file merge failed."
      end
    end
  end

  class LaneRecal
    include Pipeline::Step
    runs_tasks :count_covariates #, :plot_covariates
    runs_on :lanes

    class CountCovariates
      include Pipeline::Task
      requires_file :lane_raw_sample_bams
      dumps_file :recal_grp

      def run
	log_info "Base-quality recalibration: Count covariates"
	gatk :base_recalibrator, :knownSites => config.reference_snp_vcf, 
		:input_file => config.merged_lane_bam,
                :num_threads => nil,
		:out => config.recal_grp or error_exit "First CountCovariates failed"
      end
    end

    class PlotCovariates
      include Pipeline::Task
      requires_file :recal_grp, :merged_lane_bam
      dumps_file :recal_plot_pdf

      def run
	log_info "Generating recalibration plots."
	gatk :base_recalibrator, :knownSites => config.reference_snp_vcf, 
		:input_file => config.merged_lane_bam,
                :num_threads => nil,
                :BQSR => config.recal_grp,
		:plot_pdf_file => config.recal_plot_pdf or error_exit "Could not generate recalibration plot"
      end
    end
  end

  class LaneTableRecal
    include Pipeline::Step
    runs_tasks :table_recal
    runs_on :lanes, :chroms
    resources :threads => 1

    class TableRecal
      include Pipeline::Task
      requires_file :recal_grp, :merged_lane_bam
      dumps_file :recal_bam

      def run
	log_info "Base-quality recalibration: Table Recalibration"
	gatk :print_reads,
		:BQSR => config.recal_grp,
                :intervals => config.chrom.chrom_name,
		:input_file => config.merged_lane_bam,
                :num_threads => nil,
		:out => config.recal_bam or error_exit "TableRecalibration failed"
      end
    end
  end

  class LaneSplit
    include Pipeline::Step
    runs_tasks :split_lane, :move_files
    runs_on :lanes

    class SplitLane
      include Pipeline::Task
      requires_file :lane_recal_bams
      outs_file :lane_split_bams

      def run
	log_info "Split recalibrated bam into per-sample output bams"
	gatk :split_sam_file, :input_file => config.lane_recal_bams, :outputRoot => config.lane_split_bam_root or error_exit "Splitting bam files failed"
      end

    end
    class MoveFiles
      include Pipeline::Task
      requires_file :lane_split_bams
      outs_file :lane_sample_bams

      def run
        config.lane.samples.each do |s|
          FileUtils.mv config.lane_sample_split_bam(s), config.recal_sample_bam(s) or error_exit "Could not move realigned bam!"
        end
      end
    end
  end

  class PatientMerge
    include Pipeline::Step
    runs_tasks :merge_bam, :mark_duplicates
    runs_on :patients
    resources :threads => 1

    class MergeBam
      include Pipeline::Task
      requires_files :patient_raw_sample_bams
      dumps_file :merged_patient_bam

      def run
	log_info "Merging per-sample bam files for bulk recalibration"
        picard :merge_sam_files, :CREATE_INDEX => :true, :OUTPUT => config.merged_patient_bam, :INPUT => config.patient_raw_sample_bams or error_exit "bam file merge failed."
      end
    end

    class MarkDuplicates
      include Pipeline::Task
      requires_file :merged_patient_bam
      outs_file :raw_patient_bam, :duplication_metrics

      def run
	log_info "Mark duplicates"
	picard :mark_duplicates, :INPUT => config.merged_patient_bam,
          :OUTPUT => config.raw_patient_bam, :METRICS_FILE => config.duplication_metrics, 
          :CREATE_INDEX => :true,
          :REMOVE_DUPLICATES => :false or error_exit "Mark duplicates failed"
      end
    end
  end

  class PatientRealign
    include Pipeline::Step
    runs_tasks :create_intervals, :realign_indels
    runs_on :patients, :chroms
    resources :threads => 1

    class CreateIntervals
      include Pipeline::Task
      requires_file :raw_patient_bam
      dumps_file :patient_intervals

      def run
	log_info "Creating intervals for indel detection"
	gatk :realigner_target_creator,
          :known => config.reference_indel_vcf,
          :intervals => config.chrom.chrom_name,
          :num_threads => nil,
          :input_file => config.raw_patient_bam,
          :out => config.patient_intervals or error_exit "Interval creation failed"
      end
    end

    class RealignIndels
      include Pipeline::Task
      requires_file :raw_patient_bam, :patient_intervals
      dumps_file :realigned_patient_bam

      def run
	log_info "Indel realignment"
	gatk :indel_realigner,
			:knownAlleles => config.reference_indel_vcf,
                        :consensusDeterminationModel => :USE_READS,
                        :intervals => config.chrom.chrom_name,
			:input_file => config.raw_patient_bam,
                        :num_threads => nil,
			:targetIntervals => config.patient_intervals,
			:out => config.realigned_patient_bam or error_exit "Indel realignment failed"
      end
    end
  end

  class PatientSplit
    include Pipeline::Step
    runs_tasks :split_bam, :move_files
    runs_on :patients
    resources :threads => 1

    class SplitBam
      include Pipeline::Task
      requires_file :realigned_patient_bams
      dumps_file :patient_split_bams

      def run
	log_info "Split recalibrated bam into per-sample output bams"
	gatk :split_sam_file, :input_file => config.realigned_patient_bams, :outputRoot => config.patient_split_bam_root or error_exit "Splitting bam files failed"
      end
    end

    class MoveFiles
      include Pipeline::Task
      requires_files :patient_split_bams
      outs_file :patient_sample_bams

      def run
        config.patient.samples.each do |s|
          FileUtils.mv config.patient_sample_split_bam(s), config.sample_bam(s) or error_exit "Could not move realigned bam!"
          sam_index config.sample_bam(s) or error_exit "Indexing failed"
        end
      end
    end
  end

  class MakeSamples
    include Pipeline::Step
    runs_tasks :sample_index
    runs_on :samples

    class SampleIndex
      include Pipeline::Task
      requires_file :sample_bam
      outs_file :sample_bai

      def run
        sam_index config.sample_bam or error_exit "Indexing failed"
      end
    end
  end
end
