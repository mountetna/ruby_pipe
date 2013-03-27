#!/usr/bin/env ruby
#
def lib_file file
  File.expand_path(File.join( File.dirname(__FILE__), '../' + file ) )
end

$: << lib_file("lib")

require 'vcf'
require 'colored'
require 'net/http/persistent'
require 'resolv-replace'
require 'json'

if ARGV.length < 2
  puts "filter_pindel <config_file> <normal_name> <tumor_name> <pindel vcf>"
  exit
end

config_file = ARGV.shift
normal = ARGV.shift
tumor = ARGV.shift

puts [ :gene, :chrom, :pos, :ref_allele, :alt_allele, :tumor_ref_count, :tumor_alt_count, :normal_ref_count, :normal_alt_count, :variant_classification, :protein_change, :class ].join("\t")

ARGV.each do |f|
  vcf = VCF.new f, config_file
  vcf.each do |l|
    next if l.skip_genotype?([:pindel, :normal] => normal) || l.skip_genotype?([:pindel, :tumor] => tumor)
    next if l.skip_oncotator?

    puts [ l.onco.txp_gene, l.chrom, l.pos, l.ref, l.alt, 
      "-", l.genotype(tumor).ref_count, "-", l.genotype(normal).ref_count, 
      l.onco.txp_variant_classification, l.onco.txp_protein_change, l.onco.pph2_class ].join("\t")
  end
end
