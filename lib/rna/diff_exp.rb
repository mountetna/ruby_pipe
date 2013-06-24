#!/usr/bin/env ruby
module Rna
  class DiffExp
    include Pipeline::Step
    runs_tasks :cuff_diff_compare, :diff_exp_table
    runs_on :diff_exps
    resources :threads => 12

    class CuffDiffCompare
      include Pipeline::Task
      requires_files :replicate_bams, :normal_bams
      dumps_file :gene_exp_diff, :cuffdiff_dir

      def run
        log_info "Running cuffdiff on samples"
        cuffdiff :out => config.cuffdiff_dir, :sample1 => config.normal_bams, :sample2 => config.replicate_bams or error_exit "Could not run cuffdiff"
      end
    end
    class DiffExpTable
      include Pipeline::Task
      requires_file :gene_exp_diff
      outs_file :diff_exp_table

      def run
        diff = HashTable.new config.gene_exp_diff
        diff.select! do |l|
          l.q_value.to_f <= 0.001
        end
        diff.sort_by! {|l| l.q_value.to_f }
        diff.print config.diff_exp_table
      end
    end
  end
end
