#!/usr/bin/env ruby
require 'hash_table'
require 'sam'
module Rna
  class KallistoCount
    include Pipeline::Step
    runs_tasks :count_reads, :count_genes
    resources :threads => 2
    runs_on :samples, :replicates

    class CountReads
      include Pipeline::Task
      requires_file :input_fastq1s, :input_fastq2s
      outs_file :kallisto_abundance_tsv, :kallisto_genome_bam
      
      def run
        fastqs = config.input_fastq1s.zip(config.input_fastq2s).flatten
        kallisto_quant index: config.kallisto_index,
          output_dir: config.kallisto_scratch,
          bam: config.kallisto_genome_bam,
          input: fastqs or error_exit "Could not run kallisto"
      end
    end

    class CountGenes
      include Pipeline::Task

      requires_file :kallisto_abundance_tsv
      outs_file :kallisto_genes_results

      def run
        abundance = HashTable.new(
          types: {
            length: :float,
            eff_length: :float,
            est_counts: :float,
            tpm: :float
          }
        ).parse config.kallisto_abundance_tsv

        genes = {}

        gene_transcripts = HashTable.new(index: [ :transcript ]).parse(config.genes_transcripts)

        abundance.each do |txp|
          gene = gene_transcripts.index[:transcript][txp.target_id].first.gene
          genes[gene] ||= {
            gene_id: gene,
            est_counts: 0,
            tpm: 0
          }
          genes[gene][:est_counts] += txp.est_counts
          genes[gene][:tpm] += txp.tpm
        end

        File.open(config.kallisto_genes_results,"w") do |f|
          f.puts [ :gene_id, :est_counts, :tpm ].join("\t")
          genes.each do |gene, values|
            f.puts [ values[:gene_id], values[:est_counts], values[:tpm] ].join("\t")
          end
        end
      end
    end
  end
  class AssembleKallistoTranscripts
    include Pipeline::Step
    runs_tasks :make_coverage_table
    class MakeCoverageTable
      include Pipeline::Task
      requires_files :samples__replicates__kallisto_genes_resultss
      outs_file :kallisto_coverage_table

      def run
        names = config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r).to_sym} }.flatten
        summary = HashTable.new columns: [ :gene_id ] + names, index: [ :gene_id ]
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            gene_exps = HashTable.new.parse config.kallisto_genes_results(rep)
            name = config.sample_replicate_name(rep).to_sym
            gene_exps.each do |exp|
              if summary.index[:gene_id][exp.gene_id].count == 0
                summary << { 
                  gene_id: exp.gene_id, 
                  name => exp.est_counts
                }
              else
                entry = summary.index[:gene_id][exp.gene_id].first
                entry.update name => exp.est_counts
              end
            end
          end
        end
        summary.print config.kallisto_coverage_table
      end
    end
  end
end
