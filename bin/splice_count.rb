#!/usr/bin/env ruby
#
require 'sam'
require 'gtf'

if ARGV.length < 2
  puts "Usage: splice_count.rb <splice_sam> <splice_gtf>"
  puts "  <splice_sam> contains a list of SAM records to count"
  puts "  <splice_gtf> contains GTF lines for the gene of interest"
  exit
end

splice_sam = ARGV[0]
splice_gtf = ARGV[1]

def add_splice exon_1, exon_2
  $splices[exon_1] ||= {}
  $splices[exon_1][exon_2] ||= 0
  $splices[exon_1][exon_2] += 1
end

def count_splice exon_1, exon_2
  $splices ||= {}
  add_splice exon_1, exon_2
  add_splice exon_2, exon_1
end

def get_splice_exons splice
  exons = splice.map{|l| l.feature == "exon" ? [ l.start.to_i, l.end.to_i ] : nil }.compact.uniq.sort_by(&:first)
  exons.reverse! if splice.first.strand == "-"
  return exons
end

def count_splices sam, exons
  sam.each do |read|
    read.junctions.each do |start,stop|
      e1 = exons.index{|e| e.last == start}
      e2 = exons.index{|e| e.first == stop}
      next if !e1 || !e2
      count_splice "exon#{e1+1}", "exon#{e2+1}"
    end
  end
end

splice = GTF.new splice_gtf
sam = Sam.read splice_sam
exons = get_splice_exons splice

count_splices sam, exons

puts exons.size.times.map{|j| "exon#{j+1}"}.join("\t")
exons.size.times.each do |i|
  exon_1 = "exon#{i+1}"
  puts exon_1 + "\t" + exons.size.times.map{|j|
    exon_2 = "exon#{j+1}"
    ($splices[exon_1] ? $splices[exon_1][exon_2] : nil) || 0
  }.join("\t")
end
