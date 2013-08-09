#!/usr/bin/env ruby
require 'hash_table'
module Rna
  class DetectFusions
    include Pipeline::Step
    runs_tasks :create_input_fastq, :chimera_scan
    audit_report :sample_replicate_name
    resources :threads => 12
    runs_on :replicates

    class CreateInputFastq
      include Pipeline::Task
      requires_files :input_fastq1s, :input_fastq2s
      dumps_file :chimera_fq1, :chimera_fq2

      def run
        config.diff_exps
        run_cmd "zcat #{config.input_fastq1s.join(" ")} > #{config.chimera_fq1}" or error_exit "zcat failed"
        run_cmd "zcat #{config.input_fastq2s.join(" ")} > #{config.chimera_fq2}" or error_exit "zcat failed"
      end
    end

    class ChimeraScan
      include Pipeline::Task
      requires_files :chimera_fq1, :chimera_fq2
      outs_file :chimera_bedpe

      def run
        log_info "Running chimerascan"
        chimerascan config.chimera_fq1, config.chimera_fq2, config.chimera_dir or error_exit "Could not run chimerascan"
      end
    end
  end
end
