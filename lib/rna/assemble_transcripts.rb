#!/usr/bin/env ruby
require 'tempfile'

module Rna
  class AssembleTranscripts
    include Pipeline::Step
    runs_tasks :make_table

    class MakeTable
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
  end
end
