#!/usr/bin/env ruby
module Genome
  class LaneRecal
    include Pipeline::Step
    runs_tasks :count_covariates
    runs_on :lanes
    audit_report :lane_name

    class CountCovariates
      include Pipeline::Task
      requires_file :lane_aligned_bams
      dumps_file :recal_grp

      def run
	log_info "Base-quality recalibration: Count covariates"
	gatk :base_recalibrator, :knownSites => config.reference_snp_vcf, 
		:input_file => config.lane_aligned_bams,
                :num_threads => nil,
		:out => config.recal_grp or error_exit "First CountCovariates failed"
      end
    end
  end

  class TableRecal
    include Pipeline::Step
    runs_tasks :table_recal
    runs_on :samples, :inputs
    resources :threads => 1

    class TableRecal
      include Pipeline::Task
      requires_file :recal_grp, :input_aligned_bams
      dumps_file :recal_bam

      def run
	log_info "Base-quality recalibration: Table Recalibration"
	gatk :print_reads,
		:BQSR => config.recal_grp,
		:input_file => config.input_aligned_bams,
                :num_threads => nil,
		:out => config.recal_bam or error_exit "TableRecalibration failed"
      end
    end
  end
end
