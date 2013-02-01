#!/usr/bin/env ruby
require 'hash_table'
require 'vcf'
require 'fileutils'
require 'net/http/persistent'
require 'resolv-replace'
require 'json'

module Rna
  class FilterMuts
    include Pipeline::Step
    runs_tasks :oncotate_muts

    class OncotateMuts
      include Pipeline::Task
      requires_file :ug_filtered_vcf
      outs_file :sample_mutations

      def run 
        log_info "Filter mutations by frequency"
        vcf, normal, tumor = VCF::Sample.new(config.ug_filtered_vcf, config.mutations_config), config.normal_name, config.sample_name
        log_info "VCF loaded"

        File.open(config.sample_mutations, "w") do |f|
          f.puts [ :gene, :chrom, :pos, :ref_allele, :alt_allele, :tumor_ref_count, :tumor_alt_count, :normal_ref_count, :normal_alt_count, :variant_classification, :protein_change, :class ].join("\t")
          vcf.each do |l|
            next if l.skip_genotype?(:normal => normal) || l.skip_genotype?(:tumor => tumor)
            next if l.skip_oncotator?

            f.puts [ l.onco.txp_gene, l.chrom, l.pos, l.ref, l.alt, 
              l.genotype(tumor).ref_count, l.genotype(tumor).alt_count, l.genotype(normal).ref_count, l.genotype(normal).alt_count, 
              l.onco.txp_variant_classification, l.onco.txp_protein_change, 
              l.onco.pph2_class ].join("\t")
          end
          http.shutdown
        end
      end
    end
  end
end
