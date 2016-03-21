require 'hash_table'
require 'fileutils'
require 'vcf'
require 'maf'

module Exome
  class GermlineMutDet
    include Pipeline::Step
    runs_tasks :freebayes
    resources :threads => 1
    runs_on :normal_samples, :chroms

    class Freebayes
      include Pipeline::Task
      requires_files :sample_bam
      dumps_file :germline_unfiltered_vcf
      
      def run
        freebayes :region => config.chrom_name,
          :bam => config.sample_bam,
          :output => config.germline_unfiltered_vcf or error_exit "Could not run FreeBayes"
      end
    end
  end
  class GermlineMutFilter
    include Pipeline::Step
    runs_tasks :snp_eff_annotate_germline_vcf

    runs_on :normal_samples, :chroms
    resources memory: "4gb"

    class SnpEffAnnotateGermlineVcf
      include Pipeline::Task

      requires_file :germline_unfiltered_vcf
      dumps_file :germline_annotated_vcf

      def run
        snpeff config.germline_unfiltered_vcf, config.germline_annotated_vcf or error_exit "Could not run snpeff"
      end
    end
  end
end
