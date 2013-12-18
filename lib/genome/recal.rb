#!/usr/bin/env ruby
module Genome
 class LibraryMerge
    include Pipeline::Step
    runs_tasks :merge_bam
    resources :threads => 12, :walltime => 50

    class MergeBam
      include Pipeline::Task
      requires_files :dedup_sample_bams
      dumps_file :raw_library_bam

      def run
	log_info "Merging per-sample bam files for bulk recalibration"
        picard :merge_sam_files, :CREATE_INDEX => :true, :USE_THREADING => :true,
          :OUTPUT => config.raw_library_bam, :INPUT => config.dedup_sample_bams,
          :SORT_ORDER => :coordinate or error_exit "bam file merge failed."
      end
    end
  end

  class Recal
    include Pipeline::Step
    runs_tasks :create_intervals, :realign_indels, :count_covariates, :table_recal, :recal_index
    runs_on :chroms
    resources :threads => 1, :walltime => 50

    class CreateIntervals
      include Pipeline::Task
      requires_file :raw_library_bam
      dumps_file :merged_intervals

      def run
	log_info "Creating intervals for indel detection"
	gatk :realigner_target_creator,
          :known => config.reference_indel_vcf,
          :intervals => config.chrom,
          :num_threads => nil,
          :input_file => config.raw_library_bam,
          :out => config.merged_intervals or error_exit "Interval creation failed"
      end
    end

    class RealignIndels
      include Pipeline::Task
      requires_file :raw_library_bam, :merged_intervals
      dumps_file :realigned_bam

      def run
	log_info "Indel realignment"
	gatk :indel_realigner,
			:knownAlleles => config.reference_indel_vcf,
                        :consensusDeterminationModel => :USE_READS,
                        :intervals => config.chrom,
			:input_file => config.raw_library_bam,
                        :num_threads => nil,
			:targetIntervals => config.merged_intervals,
			:out => config.realigned_bam or error_exit "Indel realignment failed"
      end
    end

    class CountCovariates
      include Pipeline::Task
      requires_file :realigned_bam
      dumps_file :recal_grp

      def run
	log_info "Base-quality recalibration: Count covariates"
	gatk :base_recalibrator, :knownSites => config.reference_snp_vcf, 
		:input_file => config.realigned_bam,
                :num_threads => nil,
		:out => config.recal_grp or error_exit "First CountCovariates failed"
      end
    end

    class TableRecal
      include Pipeline::Task
      requires_file :recal_grp, :realigned_bam
      dumps_file :recal_bam

      def run
	log_info "Base-quality recalibration: Table Recalibration"
	gatk :print_reads,
		:BQSR => config.recal_grp,
		:input_file => config.realigned_bam,
                :num_threads => nil,
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
  end
    
  class LibrarySplit
    include Pipeline::Step
    runs_tasks :split_bam 
    resources :threads => 1, :walltime => 50

    class SplitBam
      include Pipeline::Task
      requires_file :recal_bams
      dumps_file :split_bams

      def run
	log_info "Split recalibrated bam into per-sample output bams"
	gatk :split_sam_file, :input_file => config.recal_bams, :outputRoot => config.split_bam_root or error_exit "Splitting bam files failed"
      end
    end
  end

  class MakeSamples
    include Pipeline::Step
    runs_task :sort_sample
    resources :threads => 12, :walltime => 50
    runs_on :samples

    class SortSample
      include Pipeline::Task
      requires_file :split_bam
      outs_file :sample_bam
      
      def run
	#make_samples needs more memory (-m) default 500000000
	samtools "sort -m 1000000000", config.split_bam, config.sample_bam.sub(/.bam$/,"") or error_exit "Sorting Failed"
        #sam_sort config.split_bam, config.sample_bam.sub(/.bam$/,"") or error_exit "Sorting failed"
        sam_index config.sample_bam or error_exit "Indexing failed"
      end
    end

  end
end
