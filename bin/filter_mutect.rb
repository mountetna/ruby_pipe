#!/usr/bin/env ruby
#
def lib_file file
  File.expand_path(File.join( File.dirname(__FILE__), '../' + file ) )
end

$: << lib_file("lib")

require 'mutect'
require 'colored'
require 'net/http/persistent'
require 'resolv-replace'
require 'json'


if ARGV.length < 2
  puts "filter_mutect <config_file> <mutect raw snv files>"
  exit
end

config_file = ARGV.shift
puts [ :gene, :chrom, :pos, :ref_allele, :alt_allele, :tumor_ref_count, :tumor_alt_count, :normal_ref_count, :normal_alt_count, :variant_classification, :protein_change, :class ].join("\t")
ARGV.each do |f|
  m = MuTect.new f, config_file
  m.each do |l|
    next if l.skip_mutect?
    next if l.skip_oncotator?
    puts [ l.onco.txp_gene, l.contig, l.position, l.ref_allele, l.alt_allele, 
      l.t_ref_count, l.t_alt_count, l.n_ref_count, l.n_alt_count, 
      l.onco.txp_variant_classification, l.onco.txp_protein_change, 
      l.onco.pph2_class ].join("\t")
  end
end
