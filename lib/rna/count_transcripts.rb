#!/usr/bin/env ruby
module Rna
  class CountTranscripts
    include Pipeline::Step
    runs_tasks :cufflink #, :format_transcript
    resources :threads => 12

    class Cufflink
      include Pipeline::Task
      requires_files :output_bam
      outs_file :transcripts_gtf

      def run
        log_info "Cufflinking"
        cufflinks :bam => config.output_bam, :out => config.cufflinks_scratch or error_exit "Cufflinks failed."
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
  end
end

