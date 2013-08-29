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
        combined = {}
        config.replicates.each do |rep|
          genes = HashTable.new config.gene_tracking(rep)
          genes.each do |l|
            combined[l[:gene_id]] ||= {}
            combined[l[:gene_id]][config.sample_replicate_name(rep)] = l[:FPKM]
          end
        end
        File.open config.fpkm_table, "w" do |f|
          f.puts "gene_id\t#{config.replicates.map{|r| config.sample_replicate_name(r) }.join("\t")}"
          combined.each do |gid,g|
            f.puts "#{gid}\t#{config.replicates.map{|r| g[config.sample_replicate_name(r)] }.join("\t")}"
          end
        end
      end
    end

    class MakeCoverageTable
      include Pipeline::Task
      requires_files :transcripts_covs
      outs_file :coverage_table

      def run
        combined = {}
        config.replicates.each do |rep|
          genes = HashTable.new config.transcripts_cov(rep), :header => [ :gene_id, :coverage ]
          skip = [ "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique" ]
          genes.each do |l|
            next if skip.include? l.gene_id
            combined[l.gene_id] ||= {}
            combined[l.gene_id][config.sample_replicate_name(rep)] = l.coverage
          end
        end
        File.open config.coverage_table, "w" do |f|
          f.puts "gene_id\t#{config.replicates.map{|r| config.sample_replicate_name(r) }.join("\t")}"
          combined.each do |gid,g|
            f.puts "#{gid}\t#{config.replicates.map{|r| g[config.sample_replicate_name(r)] }.join("\t")}"
          end
        end
      end
    end
  end
end
