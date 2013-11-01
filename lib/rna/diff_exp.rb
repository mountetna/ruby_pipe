#!/usr/bin/env ruby
module Rna
  class CuffDiffExp
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
          l.q_value.to_f <= config.q_value_cutoff
        end
        diff.sort_by! {|l| l.q_value.to_f }
        diff.print config.diff_exp_table
      end
    end
  end
  class DeseqDiffExp
    include Pipeline::Step
    runs_task :deseq
    runs_on :diff_exps

    class Deseq
      include Pipeline::Task
      requires_file :coverage_table
      outs_file :diff_exp_table
      def run
        r_script :deseq, :doDeseq, config.coverage_table, config.sample_name, config.normal_name, config.fdr_cutoff, config.diff_exp_table or error_exit "Could not run DESeq"
      end
    end
  end
end
