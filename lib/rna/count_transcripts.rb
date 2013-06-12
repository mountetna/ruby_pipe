#!/usr/bin/env ruby
require 'hash_table'
module Rna
  class CountTranscripts
    include Pipeline::Step
    runs_tasks :cufflink, :count_coverage #, :format_transcript
    has_tasks :cufflink, :count_coverage, :rsem_count
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
      outs_file :transcripts_gtf

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
        gene_alias = HashTable.new(config.hg19_ucsc_gene_names, :key => [ :"#kgID", :spID ])
        File.open(config.output_gtf,"w") do |f|
          File.foreach(config.transcripts_gtf) do |t|
            t.chomp!
            t.match(/gene_id "(.*?)"/) do |m|
              t.sub!(/gene_id ".*?"/,"gene_id \"#{gene_alias[m[1]][:geneSymbol]}\"") if gene_alias[m[1]]
            end
            f.puts t
          end
        end
      end
    end

    class CountCoverage
      include Pipeline::Task
      requires_file :replicate_bam
      outs_file :transcripts_cov

      def run
        log_info "Mapping coverage to reference genes"
        samtools "view -h -q 1 -F 4", config.replicate_bam, config.coverage_sam
        htseq_count :input => config.coverage_sam, :gtf => config.reference_gtf, :type => "rnaseq", :out => config.transcripts_cov or error_exit "Computing normal coverage failed."
        File.unlink config.coverage_sam
        #coverage_bed config.sample_bam, config.reference_gtf, config.normal_cov or error_exit "Computing normal coverage failed."
      end
    end
  end
end
