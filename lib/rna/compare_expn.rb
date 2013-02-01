#!/usr/bin/env ruby
module Rna
  class CompareExpn
    include Pipeline::Step
    runs_tasks :cufflink

    class Cufflink
      include Pipeline::Task
      requires_files :sorted_bam
      outs_file :transcripts_gtf

      def run
        log_info "Cufflinking"
        cufflinks config.output_bam, config.cufflinks_output or error_exit "Cufflinks failed."
      end
    end
  end
end

