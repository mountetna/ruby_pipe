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

vcf,normal,tumor = VCF::Sample.new(ARGV[0], lib_file('config/rna_mutations.yml') ),ARGV[1],ARGV[2]

puts "VCF loaded".green

vcf.each do |l|
  next if l.skip_genotype?(:normal => normal) || l.skip_genotype?(:tumor => tumor)
  next if l.skip_oncotator?

  puts [ l.onco.txp_gene, l.chrom, l.pos, l.ref, l.alt, 
    l.genotype(tumor).ref_count, l.genotype(tumor).alt_count, l.genotype(normal).ref_count, l.genotype(normal).alt_count, 
    l.onco.txp_variant_classification, l.onco.txp_protein_change, l.onco.pph2_class ].join("\t")
  #puts [ l.genotype(normal).gt, l.genotype(tumor).gt, l.chrom, l.pos, l.ref, l.alt, l.genotype(tumor).ref_count, l.genotype(tumor).alt_count, l.genotype(normal).ref_count, l.genotype(normal).alt_count, l.genotype(tumor).alt_freq, l.genotype(normal).alt_freq ].join("\t")
end
