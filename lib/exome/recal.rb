#!/usr/bin/env ruby
module Exome
  class Recal
    include Pipeline::Step
    runs_tasks :merge_bam, :create_intervals, :realign_indels, :fix_mates, :mark_duplicates, :count_covariates, :table_recal, :recal_index, :split_bam 
    class MergeBam
      include Pipeline::Task
      requires_files :recal_bams
      dumps_file :merged_bam

      def run
	log_info "Merging bam files"
        picard :merge_sam_files, :CREATE_INDEX => :true, :OUTPUT => config.merged_bam, :INPUT => config.recal_bams, :SORT_ORDER => :coordinate or error_exit "bam file merge failed."
	sam_index config.merged_bam or error_exit "First indexing failed"
      end
    end
    class CreateIntervals
      include Pipeline::Task
      requires_file :merged_bam
      dumps_file :merged_intervals

      def run
	log_info "Creating intervals for indel detection"
	gatk :realigner_target_creator,
          :known => config.thousand_genomes_indels,
          :num_threads => 24,
          :input_file => config.merged_bam,
          :out => config.merged_intervals or error_exit "Interval creation failed"
      end
    end

    class RealignIndels
      include Pipeline::Task
      requires_file :merged_bam, :merged_intervals
      dumps_file :realigned_bam

      def run
	log_info "Indel realignment"
	gatk :indel_realigner,
			:knownAlleles => config.thousand_genomes_indels,
                        :consensusDeterminationModel => :USE_READS,
			:input_file => config.merged_bam,
			:targetIntervals => config.merged_intervals,
			:out => config.realigned_bam or error_exit "Indel realignment failed"
      end
    end

    class FixMates
      include Pipeline::Task
      requires_file :realigned_bam
      dumps_file :mated_bam

      def run
	log_info "Fix mate information"
	picard :fix_mate_information, :INPUT => config.realigned_bam, :OUTPUT => config.mated_bam, :SO => :coordinate or error_exit "Fix mate information failed"
      end
    end

    class MarkDuplicates
      include Pipeline::Task
      requires_file :mated_bam
      dumps_file :dedup_bam

      def run
	log_info "Mark duplicates"
	picard :mark_duplicates, :INPUT => config.mated_bam, :OUTPUT => config.dedup_bam, :METRICS_FILE => config.recal_metrics, :REMOVE_DUPLICATES => :true, :CREATE_INDEX => :true or error_exit "Mark duplicates failed"
      end
    end

    class CountCovariates
      include Pipeline::Task
      requires_file :dedup_bam
      dumps_file :recal_grp

      def run
	log_info "Base-quality recalibration: Count covariates"
	gatk :base_recalibrator, :knownSites => config.dbsnp_vcf, 
		:input_file => config.dedup_bam,
		:out => config.recal_grp or error_exit "First CountCovariates failed"
      end
    end

    class TableRecal
      include Pipeline::Task
      requires_file :recal_grp
      dumps_file :recal_bam

      def run
	log_info "Base-quality recalibration: Table Recalibration"
	gatk :print_reads,
		:BQSR => config.recal_grp,
		:input_file => config.dedup_bam,
		:out => config.recal_bam or error_exit "TableRecalibration failed"
      end
    end

    class RecalIndex
      include Pipeline::Task
      requires_file :recal_bam
      dumps_file :recal_bai

      def run
	log_info "Index recal bam file"
	sam_index config.recal_bam or error_exit "Recal indexing failed"
      end
    end
    
    class SplitBam
      include Pipeline::Task
      requires_file :recal_bam, :recal_bai
      dumps_file :split_bams
      outs_file :output_bams

      def run
	log_info "Split recalibrated bam into per-sample output bams"
	gatk :split_sam_file,
		:input_file => config.recal_bam,
		:outputRoot => config.split_bam or error_exit "Splitting bam files failed"

	# split off the list of files based on the BASE list, to verify
        config.sample_names.each do |sample|
          tmp = config.split_bam(sample)
          error_exit "No temp bam for #{sample} found" if !File.exists?(tmp)
	  log_info "Processing #{sample}"
          sam_sort tmp, config.output_bam(sample).sub(/.bam$/,"") or error_exit "Sorting #{tmp} failed"
          sam_index config.output_bam(sample) or error_exit "Sorting #{sample} output failed"
        end
      end
    end
  end
end
