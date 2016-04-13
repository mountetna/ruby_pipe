#!/usr/bin/env ruby
require 'tempfile'

module Rna
  class AssembleTranscripts
    include Pipeline::Step
    runs_tasks :make_fpkm_table, :make_coverage_table

    class MakeFpkmTable
      include Pipeline::Task
      requires_files :gene_trackings
      outs_file :fpkm_table

      def run
        summary = HashTable.new columns: [ :gene_id ] + config.samples__replicatess.map{ |rep| config.sample_replicate_name(rep).to_sym }, index: [ :gene_id ]
        config.samples__replicatess.each do |rep|
          gene_exps = HashTable.new.parse config.gene_tracking(rep)
          name = config.sample_replicate_name(rep).to_sym
          gene_exps.each do |exp|
            if summary.index[:gene_id][exp.gid].count == 0
              summary << { 
                gene_id: exp.gene_id, 
                name => gene.fpkm
              }
            else
              entry = summary.index[:gene_id][exp.gene_id].first
              entry.update name => exp.fpkm
            end
          end
        end
        summary.print config.fpkm_table
      end
    end

    class MakeCoverageTable
      include Pipeline::Task
      requires_files :transcripts_covs
      outs_file :coverage_table

      def run
        combined = {}
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            genes = HashTable.new config.transcripts_cov(rep), :header => [ :gene_id, :coverage ]
            skip = [ "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique" ]
            genes.each do |l|
              next if skip.include? l.gene_id
              combined[l.gene_id] ||= {}
              combined[l.gene_id][config.sample_replicate_name(rep)] = l.coverage
            end
          end
        end
        File.open config.coverage_table, "w" do |f|
          f.puts "gene_id\t#{config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r) } }.flatten.join("\t")}"
          combined.each do |gid,g|
            f.puts "#{gid}\t#{config.samples.map{|s| s.replicates.map{|r| g[config.sample_replicate_name(r)] } }.flatten.join("\t")}"
          end
        end
      end
    end
  end
  class AssembleRsemTranscripts
    include Pipeline::Step
    runs_tasks :make_tpm_table, :make_coverage_table

    class MakeTpmTable
      include Pipeline::Task
      requires_files :samples__replicates__rsem_genes_resultss
      outs_file :tpm_table

      def run
        names = config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r).to_sym} }.flatten
        summary = HashTable.new columns: [ :gene_id ] + names, index: [ :gene_id ]
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            gene_exps = HashTable.new.parse config.rsem_genes_results(rep)
            name = config.sample_replicate_name(rep).to_sym
            gene_exps.each do |exp|
              if summary.index[:gene_id][exp.gene_id].count == 0
                summary << { 
                  gene_id: exp.gene_id, 
                  name => exp.tpm
                }
              else
                entry = summary.index[:gene_id][exp.gene_id].first
                entry.update name => exp.tpm
              end
            end
          end
        end
        summary.print config.tpm_table
      end
    end

    class MakeCoverageTable
      include Pipeline::Task
      requires_files :samples__replicates__rsem_genes_resultss
      outs_file :coverage_table

      def run
        names = config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r).to_sym} }.flatten
        summary = HashTable.new columns: [ :gene_id ] + names, index: [ :gene_id ]
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            gene_exps = HashTable.new.parse config.rsem_genes_results(rep)
            name = config.sample_replicate_name(rep).to_sym
            gene_exps.each do |exp|
              if summary.index[:gene_id][exp.gene_id].count == 0
                summary << { 
                  gene_id: exp.gene_id, 
                  name => exp.expected_count
                }
              else
                entry = summary.index[:gene_id][exp.gene_id].first
                entry.update name => exp.expected_count
              end
            end
          end
        end
        summary.print config.coverage_table
      end
    end
  end
end
