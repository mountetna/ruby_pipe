#!/usr/bin/env ruby
module Exome
  class CopyNumber
    include Pipeline::Step
    runs_tasks :compute_coverage
    class ComputeCoverage
      include Pipeline::Task
      requires_files :tumor_bam, :normal_bam, :interval_bed
      dumps_files :tumor_cov, :normal_cov

      def run
        if !File.exists? config.normal_cov
          coverage_bed config.normal_bam, config.interval_bed, config.normal_cov or error_exit "Computing normal coverage failed."
        end
        coverage_bed config.tumor_bam, config.interval_bed, config.tumor_cov or error_exit "Computing tumor coverage failed."
      end
    end
  end
end
