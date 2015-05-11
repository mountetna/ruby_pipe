#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'vcf'

module Exome
  class RunAscat
    include Pipeline::Step
    runs_on :tumor_samples
    audit_report :sample_name, :normal_name
    runs_tasks :compute_baf, :ascat_purity_ploidy

    class ComputeBaf
      include Pipeline::Task
      requires_file :ug_snps_vcf
      outs_file :tumor_baf, :normal_baf

      def run
        vcf = VCF.read config.ug_snps_vcf
        tb = HashTable.new :columns => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        nb = HashTable.new :columns => [ :chromosome, :position, :BAF, :Alt_count, :Tot_count ]
        vcf.each do |v|
          next if !v.genotype(config.normal.sample_id).heterozygous?
          next if !v.genotype(config.sample.sample_id).callable?
          next if v.filter != "PASS"
          next if v.id[0..1] != "rs"
          tb.add_line :chromosome => v.seqname, :position => v.pos,
            :BAF => v.genotype(config.sample.sample_id).alt_freq,
            :Alt_count => v.genotype(config.sample.sample_id).alt_count,
            :Tot_count => v.genotype(config.sample.sample_id).depth
          nb.add_line :chromosome => v.seqname, :position => v.pos,
            :BAF => v.genotype(config.normal.sample_id).alt_freq,
            :Alt_count => v.genotype(config.normal.sample_id).alt_count,
            :Tot_count => v.genotype(config.normal.sample_id).depth
        end
        tb.print config.tumor_baf
        nb.print config.normal_baf
      end
    end

    class AscatPurityPloidy
      include Pipeline::Task
      requires_file :sample_exon_cnr, :tumor_baf, :normal_exon_cnr, :normal_baf, :interval_bed
      outs_file :tumor_ascat_rdata, :tumor_ascat_txt

      def run
        r_script :segment, :doAscatPurityPloidy, config.tumor_baf, config.normal_baf, config.sample_exon_cnr, config.normal_exon_cnr, config.interval_bed, config.tumor_ascat_rdata, config.tumor_ascat_txt or error_exit "ASCAT failed"
      end
    end
  end

  class RunAbsolute
    include Pipeline::Step
    runs_on :tumor_samples
    audit_report :sample_name, :normal_name
    runs_tasks :absolute_purity_ploidy
    class AbsolutePurityPloidy
      include Pipeline::Task
      requires_file :sample_exon_cnr, :all_muts_maf
      outs_file :absolute_rdata, :absolute_pdf

      def run
        r_script :absolute, :callSample, config.sample_name, config.tumor_cnr_seg, config.all_muts_maf, config.absolute_scratch or error_exit "Absolute failed"

        FileUtils.cp config.absolute_scratch_pdf, config.absolute_pdf if File.exists? config.absolute_scratch_pdf
        FileUtils.cp config.absolute_scratch_rdata, config.absolute_rdata if File.exists? config.absolute_scratch_rdata
      end
    end
  end

  class ReviewAbsolute
    include Pipeline::Step
    runs_tasks :create_review_object, :extract_review_results

    class CreateReviewObject
      include Pipeline::Task

      requires_files :absolute_rdatas
      dumps_file :review_table

      def run
        r_script :absolute, :createReview, config.cohort_name, config.absolute_review_dir, config.absolute_rdatas.join(",") or error_exit "Absolute failed"
      end
    end
    class ExtractReviewResults
      include Pipeline::Task

      requires_file :reviewed_table, :absolute_modes
      outs_file  :absolute_calls

      def run
        r_script :absolute, :extractReview, config.reviewed_table, "Exome.Pipeline", config.absolute_modes, config.absolute_review_dir, config.cohort_name or error_exit "Absolute failed"
        FileUtils.mv config.absolute_calls_scratch, config.absolute_calls 
        config.tumor_samples.each do |sample|
          FileUtils.mv config.absolute_segs_scratch(sample), config.tumor_absolute_seg(sample)
        end
      end
    end
  end
end
