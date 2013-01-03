#!/usr/bin/env ruby
module Exome
  class CopyNumber
    include Pipeline::Step
    runs_tasks :compute_coverage, :compute_ratio
    class ComputeCoverage
      include Pipeline::Task
      requires_files :tumor_bam, :normal_bam, :interval_list
      dumps_files :tumor_cov, :normal_cov

      def run
        coverage_bed config.tumor_bam, config.interval_list, config.tumor_cov or error_exit "Computing tumor coverage failed."
        coverage_bed config.normal_bam, config.interval_list, config.normal_cov or error_exit "Computing normal coverage failed."
      end
    end
  end
end
