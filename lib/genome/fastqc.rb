#!/usr/bin/env ruby
module Genome
  class FastQc
    include Pipeline::Step
    runs_tasks :run_fast_qc, :make_fastqc_pdf, :verify_fastqc
    runs_on :samples, :fastqs
    # resources :threads => 12

    class RunFastQc
      include Pipeline::Task
      requires_file :fastq_file
      outs_file :fastqc_summary, :fastqc_html

      def run
        fastqc config.fastq_file, config.sample_scratch or error_exit "Could not deal with fastq"
      end
    end

    class VerifyFastqc
      include Pipeline::Task
      requires_file :fastqc_summary

      def fastqc_summary_fails? file
        File.foreach(file)
          .map{|l| l.split(/\t/).first }
          .count{|result| result == "FAIL"} > 1
      end

      def run
        if config.verify_fastq_quality
          error_exit "Fastqc failed quality check!" if fastqc_summary_fails?(config.fastqc_summary)
        end
      end
    end

    class MakeFastqcPdf
      include Pipeline::Task
      requires_file :fastqc_html
      outs_file :fastqc_pdf

      def run
        htmldoc config.fastqc_html, config.fastqc_pdf or error_exit "Could not make fastqc pdf"
      end
    end
  end
end
