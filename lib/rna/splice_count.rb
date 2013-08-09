#!/usr/bin/env ruby
require 'gtf'
require 'sam'
module Rna
  class SpliceCount
    include Pipeline::Step
    runs_tasks :create_read_cigar, :count_exon_splices
    audit_report :sample_replicate_name
    runs_on :replicates

    class CreateReadCigar
      include Pipeline::Task
      requires_file :replicate_bam, :splice_gtf
      outs_file :splice_sam
      
      def run
        splice = GTF.new config.splice_gtf
        # get the max and min intervals
        min = splice.map{|s| s.start.to_i}.min
        max = splice.map{|s| s.end.to_i}.max
        chr = splice.first.seqname

        run_cmd "samtools view -q 20 #{config.replicate_bam} #{chr}:#{min}-#{max} > #{config.splice_sam}" or error_exit "Could not get reads for the gene interval"
      end
    end

    class CountExonSplices
      include Pipeline::Task
      requires_file :splice_sam, :splice_gtf
      outs_file :exon_splice_counts

      def add_splice exon_1, exon_2
        @splices[exon_1] ||= {}
        @splices[exon_1][exon_2] ||= 0
        @splices[exon_1][exon_2] += 1
      end

      def count_splice exon_1, exon_2
        @splices ||= {}
        add_splice exon_1, exon_2
        add_splice exon_2, exon_1
      end

      def run
        splice = GTF.new config.splice_gtf
        sam = Sam.read config.splice_sam

        exons = splice.map{|l| l.feature == "exon" ? [ l.start.to_i, l.end.to_i ] : nil }.compact.uniq.sort_by(&:first)
        exons.reverse! if splice.first.strand == "-"

        sam.each do |read|
          read.junctions.each do |start,stop|
            e1 = exons.index{|e| e.last == start}
            e2 = exons.index{|e| e.first == stop}
            next if !e1 || !e2
            count_splice "exon#{e1+1}", "exon#{e2+1}"
          end
        end

        File.open(config.exon_splice_counts,"w") do |f|
          f.puts exons.size.times.map{|j| "exon#{j+1}"}.join("\t")
          exons.size.times.each do |i|
            exon_1 = "exon#{i+1}"
            f.puts exon_1 + "\t" + exons.size.times.map{|j|
              exon_2 = "exon#{j+1}"
              (@splices[exon_1] ? @splices[exon_1][exon_2] : nil) || 0
            }.join("\t")
          end
        end

        #py_script :count_exonjunctions, :args => [ config.splice_sam, config.splice_gtf, strand ], :out => config.exon_splice_counts or error_exit "Exon counting failed"
      end
    end
  end
end
