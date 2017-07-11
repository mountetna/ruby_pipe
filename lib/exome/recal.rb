#!/usr/bin/env ruby
module Exome
  # the Broad folks do:
  # 1. Merge all BAMs for a lane (this might not be possible, but do your best) 
  # 2. recalibrate per lane
  # 3. Split and merge bams for each patient
  # 4. Dedup and indel realign
  # 5. Split and write bams
  class LaneRecal
    include Pipeline::Step
    runs_tasks :count_covariates
    runs_on :lanes
    resources memory: "10gb"

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
    resources memory: "10gb"

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
end
