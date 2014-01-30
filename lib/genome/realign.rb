#!/usr/bin/env ruby
module Genome
  class PatientRealign
    include Pipeline::Step
    runs_tasks :create_intervals, :realign_indels, :split_bam
    runs_on :patients, :chroms
    audit_report :patient_name
    resources :threads => 1

    class CreateIntervals
      include Pipeline::Task
      requires_file :patient_recal_bams
      dumps_file :patient_intervals

      def run
	log_info "Creating intervals for indel detection"
	gatk :realigner_target_creator,
          :known => config.reference_indel_vcf,
          :intervals => config.chrom.chrom_name,
          :num_threads => nil,
          :input_file => config.patient_recal_bams,
          :out => config.patient_intervals or error_exit "Interval creation failed"
      end
    end

    class RealignIndels
      include Pipeline::Task
      requires_file :patient_recal_bams, :patient_intervals
      dumps_file :realigned_patient_bam

      def run
	log_info "Indel realignment"
	gatk :indel_realigner,
			:knownAlleles => config.reference_indel_vcf,
                        :consensusDeterminationModel => :USE_READS,
                        :intervals => config.chrom.chrom_name,
			:input_file => config.patient_recal_bams,
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
