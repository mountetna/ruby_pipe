#!/usr/bin/env ruby
require 'hash_table'
module Rna
  class CountTranscripts
    include Pipeline::Step
    runs_tasks :cufflink, :format_transcript, :sort_sam, :count_coverage
    has_tasks :cufflink, :format_transcript, :sort_sam, :count_coverage, :rsem_count
    audit_report :sample_replicate_name
    resources :threads => 12
    runs_on :replicates

    class RsemCount
      include Pipeline::Task
      requires_file :replicate_bam
      outs_file :rsem_genes_results
      
      def run
        rsem :calculate_expression, :input => config.replicate_bam, :sample => "#{config.sample_name}.#{config.job_item.replicate_name}", :output => config.rsem_output_dir, :bam => true
      end
    end

    class Cufflink
      include Pipeline::Task
      requires_files :replicate_bam
      dumps_file :transcripts_gtf

      def run
        log_info "Running cufflinks"
        cufflinks :bam => config.replicate_bam, :out => config.cufflinks_scratch or error_exit "Cufflinks failed."
      end
    end

    class FormatTranscript
      include Pipeline::Task
      requires_file :transcripts_gtf
      outs_file :output_gtf

      def run
        log_info "Formatting transcript gtf"
        FileUtils.cp config.transcripts_gtf, config.output_gtf
      end
    end

    class SortSam
      include Pipeline::Task
      requires_file :replicate_bam
      dumps_file :coverage_sam

      def run
        run_cmd "samtools view -b -h -q 1 -f 2 #{config.replicate_bam} | samtools sort -on - #{config.coverage_sam.sub(/.sam/,"")} | samtools view -h - > #{config.coverage_sam}" or error_exit "Could not sort bam correctly."
      end
    end

    class CountCoverage
      include Pipeline::Task
      requires_file :coverage_sam
      outs_file :transcripts_cov

      def run
        log_info "Mapping coverage to reference genes"
        htseq_count :input => config.coverage_sam, :gtf => config.reference_gtf, :type => "exon", :out => config.transcripts_cov or error_exit "Computing normal coverage failed."
      end
    end
  end
end
