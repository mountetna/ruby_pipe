#!/usr/bin/env ruby
require 'maf'
require 'fasta'

module Exome
  class MutSpectrum
    include Pipeline::Step
    runs_tasks :compute_context, :plot_spectrum
    resources :threads => 1
    runs_on :tumor_samples

    class ComputeContext
      include Pipeline::Task
      requires_files :tumor_maf
      dumps_files :tumor_context_maf

      def run
	log_info "Computing context for mutations"
        m = Maf.new.parse config.tumor_maf
        m.add_column :context
        m.each do |mut|
          next if mut.variant_type != "SNP"
          ctxt = Fasta.default.get_seq mut.short_chrom, mut.start-1, mut.stop+1
          mut.context = "#{ctxt[0]}[#{mut.reference_allele}/#{mut.tumor_seq_allele1}]#{ctxt[2]}"
        end

        m.write config.tumor_context_maf
      end
    end

    class PlotSpectrum
      include Pipeline::Task

      requires_files :tumor_context_maf
      dumps_file :stratton_plot_pdf

      def run
        r_script :stratton, :doStrattonPlot, config.tumor_context_maf, config.stratton_plot_pdf
      end
    end
  end
end
