#!/usr/bin/env ruby
module Exome
  # the Broad folks do:
  # 1. Merge all BAMs for a lane (this might not be possible, but do your best) 
  # 2. recalibrate per lane
  # 3. Split and merge bams for each patient
  # 4. Dedup and indel realign
  # 5. Split and write bams
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
    runs_tasks :merge_files, :sample_index
    runs_on :samples

    class MoveFiles
      include Pipeline::Task
      requires_file :sample_split_bam
      outs_file :sample_bam

      def run
        FileUtils.mv config.sample_split_bam, config.sample_bam or error_exit "Could not move file"
      end
    end

    class MergeFiles
      include Pipeline::Task
      requires_file :sample_patient_bams
      outs_file :sample_bam

      def run
        picard :merge_sam_files, :CREATE_INDEX => :true, :OUTPUT => config.sample_bam, :INPUT => config.sample_patient_bams or error_exit "bam file merge failed."
      end
    end

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
