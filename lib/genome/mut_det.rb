#!/usr/bin/env ruby

require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

module Genome
  class VariantDet
    include Pipeline::Step
    runs_tasks :variant_caller, :variant_annot, :quality_filter
    resources :threads => 1, :walltime => 50
    runs_on :patients, :chroms

    class VariantCaller
      include Pipeline::Task
      requires_files :samples__sample_bams
      dumps_files :ug_raw_vcf
      def run
        log_info "Running GATK UnifiedGenotypeCaller on normal/tumor sample"
        gatk :unified_genotyper, 
             :genotype_likelihoods_model => :BOTH,
             :genotyping_mode => :DISCOVERY,
             :intervals => config.chrom.chrom_name,
             :input_file => config.samples__sample_bams, 
             :"dbsnp" => config.reference_snp_vcf,
             :"standard_min_confidence_threshold_for_calling" => 30.0,
             :"standard_min_confidence_threshold_for_emitting" => 10.0,
             :"min_base_quality_score" => 20,
             :"output_mode" => "EMIT_VARIANTS_ONLY",
             :"out" => config.ug_raw_vcf or error_exit "GATK UnifiedGenotypeCaller failed"
      end
    end

    class VariantAnnot
      include Pipeline::Task
      requires_file :ug_raw_vcf
      outs_files :ug_annotated_vcf 

      def run
        log_info "Running GATK VariantRecal on VCFs"
        gatk :variant_annotator, :input_file => config.normal_bam, 
             :intervals => config.chrom.chrom_name, 
             :variant => config.ug_raw_vcf, :num_threads => 1,
             :dbsnp => config.reference_snp_vcf, 
             :baq =>  :CALCULATE_AS_NECESSARY,
             :annotation => [ :QualByDepth, :RMSMappingQuality,
                              :MappingQualityZero, :LowMQ,
                              :MappingQualityRankSumTest, :FisherStrand,
                              :HaplotypeScore, :ReadPosRankSumTest, 
                              :Coverage ],
             :out => config.ug_annotated_vcf or error_exit "GATK VariantRecalibrator failed" 
      end
    end

    class QualityFilter
      include Pipeline::Task
      requires_files :ug_annotated_vcf
      dumps_file :ug_filtered_vcf

      def run
        log_info "Filtering Unified Genotyper SNPs"
        gatk :variant_filtration,
                :variant => config.ug_annotated_vcf,
                :intervals => config.chrom.chrom_name,
                :num_threads => 1,
                :baq => :CALCULATE_AS_NECESSARY,
                :filterExpression => [ '"QD < 2.0"',
                                       '"MQ < 40.0"',
                                       '"FS > 60.0"',
                                       '"HaplotypeScore > 13.0"',
                                       '"MQRankSum < -12.5"',
                                       '"ReadPosRankSum < -8.0"' ],
                :filterName => [ :QDFilter, :MQFilter, :FSFilter,
                                 :HaplotypeScoreFilter, :MQRankSumFilter,
                                 :ReadPosFilter ],
                :out => config.ug_filtered_vcf or error_exit "Unified Genotyper SNP filtration failed"

      end
    end
  end

  class MergeVariants
    include Pipeline::Step
    runs_tasks :merge_variants
    resources :threads => 12, :walltime => 50
    runs_on :patients

    class MergeVariants
      include Pipeline::Task
      requires_files :chroms__ug_filtered_vcfs
      dumps_file :ug_vcf

      def run
        log_info "Combining annotated, filatered VCFs"
        gatk :combine_variants, :variant => config.chroms__ug_filtered_vcfs, :num_threads => 1,
             :out => config.ug_vcf, :"genotypemergeoption" => "UNSORTED" or error_exit "Merging UnifiedGenotyper failed"
      end
    end
  end 

  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :pindel, :make_pindel_vcf, :patch_pindel_vcf
    has_tasks :mutect, :pindel, :make_pindel_vcf, :patch_pindel_vcf, :indelocator, :patch_indelocator_vcf, :somatic_indel_detector, :patch_somatic_indel_vcf
    resources :threads => 1, :walltime => 50
    runs_on :tumor_samples, :chroms

    def vacuum
      config.sample_names.each do |s|
        FileUtils.rm Dir.glob("#{config.scratch_dir}/#{s}/AnnoVar*")
      end
    end

    class Mutect
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_files :mutect_snv, :mutect_coverage

      def run
	log_info "Running muTect for tumor #{config.sample_name}, normal #{config.normal_name}"
        mutect "input_file:normal" => config.normal_bam, "input_file:tumor" => config.tumor_bam,
	  :intervals => config.chrom.chrom_name, 
          :out => config.mutect_snv_tmp, 
          :coverage_file => config.mutect_coverage or error_exit "muTect failed"

        FileUtils.mv config.mutect_snv_tmp, config.mutect_snv
      end
    end
    #CHANGMT: this is for BWA-aligned sequences
    class Pindel
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_files :pindel_snv_d

      def run
        log_info "Running pindel"
        pindel :bams => [
            { :bam => config.tumor_bam, :name => config.sample_name},
            { :bam => config.normal_bam, :name => config.normal_name } ],
          :tempfile => config.pindel_list,
          :chromosome => config.chrom.chrom_name, :"output-prefix" => config.pindel_snvs or error_exit "Pindel failed"
      end
    end

    class MakePindelVcf
      include Pipeline::Task
      dumps_file :pindel_unpatched_vcf

      def run
        pindel_to_vcf :pindel_output_root => config.pindel_snvs, :vcf => config.pindel_unpatched_vcf or error_exit "Pindel2VCF failed"
      end
    end

    class PatchPindelVcf
      include Pipeline::Task
      requires_file :pindel_unpatched_vcf, :normal_bam, :tumor_bam
      dumps_file :pindel_vcf

      def run
        unpatched = VCF.read config.pindel_unpatched_vcf
        unpatched.preamble_lines.find{|i| i=~ /ID=AD/}.replace "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n"
        log_info "Trying to find genotyes for #{config.normal_name}, #{config.sample_name}"
        unpatched.each do |l|
          #skip invalid ones
          if l.alt == "<INS>"
            l.invalid = true 
            next
          end
          l.format.push :DP 
          # get the actual depth for this locus
          normal_depth = count_depth config.normal_bam, l.mutation[:chrom], l.mutation[:pos]
          tumor_depth = count_depth config.tumor_bam, l.mutation[:chrom], l.mutation[:pos]
          support = l.genotype(config.normal_name).info[:AD].to_i
          if support > normal_depth
            log_error "Mpileup depth is less than pindel supporting read count!"
            log_error "#{support} #{normal_depth} #{l.chrom}:#{l.pos}" 
            normal_depth = support
          end
          l.genotype(config.normal_name).info[:DP] = normal_depth
          l.genotype(config.normal_name).info[:AD] = "#{normal_depth - support},#{support}"

          support = l.genotype(config.sample_name).info[:AD].to_i
          if support > tumor_depth
            log_error "Mpileup depth is less than pindel supporting read count!"
            log_error "#{support} #{tumor_depth} #{l.chrom}:#{l.pos}"
            tumor_depth = support
          end
          l.genotype(config.sample_name).info[:DP] = tumor_depth
          l.genotype(config.sample_name).info[:AD] = "#{tumor_depth - support},#{support}"
        end
        unpatched.write config.pindel_vcf
      end
    end

    class Indelocator
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      outs_files :indelocator_unpatched_vcf, :indelocator_metrics

      def run
        indelocator :somatic => true, :"input_file:normal" => config.normal_bam, :"input_file:tumor" => config.tumor_bam, :verboseOutput => config.indelocator_output, :metrics_file => config.indelocator_metrics, :out => config.indelocator_unpatched_vcf, :intervals => config.chrom.chrom_name or error_exit "Indelocator failed"
      end
    end
    class PatchIndelocatorVcf
      include Pipeline::Task
      requires_file :indelocator_unpatched_vcf
      dumps_file :indelocator_vcf

      def run
        unpatched = VCF.read config.indelocator_unpatched_vcf
        unpatched.preamble_lines.push "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n"
        unpatched.each do |l|
          l.format.push :DP, :AD
          l.genotype(config.normal_name).info[:DP] = l.info[:N_DP]
          l.genotype(config.normal_name).info[:AD] = [ l.info[:N_DP].to_i - l.info[:N_AC].split(/,/).first.to_i,
                                                       l.info[:N_AC].split(/,/).first.to_i ].join(",")

          l.genotype(config.sample_name).info[:DP] = l.info[:T_DP]
          l.genotype(config.sample_name).info[:AD] = [ l.info[:T_DP].to_i - l.info[:T_AC].split(/,/).first.to_i,
                                                       l.info[:T_AC].split(/,/).first.to_i ].join(",")
        end
        unpatched.write config.indelocator_vcf
      end
    end

    class SomaticIndelDetector
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      outs_files :somaticindel_unpatched_vcf, :somaticindel_verbose

      def run
        somatic_indel_detector :"input_file:normal" => config.normal_bam, 
          :"input_file:tumor" => config.tumor_bam, 
          :verboseOutput => config.somaticindel_verbose, 
          :out => config.somaticindel_unpatched_vcf, 
          :filter_expressions => "'T_COV<6||N_COV<4||T_INDEL_CF<0.7'",
          :intervals => config.chrom.chrom_name or error_exit "SomaticIndelDetector failed"
      end
    end

    class SomaticIndelOut
      include Enumerable
      FIELDS = %r{
        (?<field_name> [A-Z_]+ ){0}
        (?<params> [A-Z\/]+ ){0}
        (?<fields> [\d\.\/-]+){0}
        \A\g<field_name>(\[\g<params>\])?:?\g<fields>?\Z
      }x
      FIELD_KEYS = {
         :N_OBS_COUNTS => [ :alt, :other, :total],
         :N_AV_MM => [:alt,:ref],
         :N_AV_MAPQ => [:alt,:ref],
         :N_NQS_MM_RATE => [:alt,:ref],
         :N_NQS_AV_QUAL => [:alt,:ref],
         :N_STRAND_COUNTS => [:alt_forward,:alt_reverse,:ref_forward,:ref_reverse],
         :T_OBS_COUNTS => [:alt,:other,:total],
         :T_AV_MM => [:alt,:ref],
         :T_AV_MAPQ => [:alt,:ref],
         :T_NQS_MM_RATE => [:alt,:ref],
         :T_NQS_AV_QUAL => [:alt,:ref],
         :T_STRAND_COUNTS => [:alt_forward,:alt_reverse,:ref_forward,:ref_reverse],
      }
      def initialize file
        @muts = Hash[File.foreach(file).map do |l|
           mut = SomaticIndelOut::Indel.new l.chomp.split(/\t/)
           [ mut.key, mut ]
        end]
      end

      class Indel
        def initialize arr
          @mut = {
            :chrom => arr.shift,
            :start => arr.shift,
            :stop => arr.shift,
            :indel => arr.shift
          }
          read_fields arr
        end

        def read_fields arr
          arr.each do |blob|
            r = blob.match(FIELDS)
            name = r[:field_name].downcase.to_sym
            @mut[ r[:field_name].downcase.to_sym ] = case
            when r[:fields] && (FIELD_KEYS[name] || r[:params])
              Hash[
                (FIELD_KEYS[name] || r[:params].split(%r!/!)).zip(
                  r[:fields].split(%r!/!)
                )]
            when r[:fields] && !FIELD_KEYS[name] && !r[:params]
              r[:fields]
            else
              true
            end
          end
        end

        def key
          [ @mut[:chrom], @mut[:start].to_i, @mut[:stop].to_i ]
        end

        def method_missing meth, *args, &block
          @mut[meth] || super(meth, *args, &block)
        end
        
        def normal_depth
          @normal_depth ||= n_obs_counts[:total].to_i
        end
        def normal_alt_count
          @normal_alt ||= n_obs_counts[:alt].to_i
        end
        def normal_allelic_depth
          [ normal_depth - normal_alt_count, normal_alt_count ].join ","
        end

        def tumor_depth
          @tumor_depth ||= t_obs_counts[:total].to_i
        end
        def tumor_alt_count
          @tumor_alt ||= t_obs_counts[:alt].to_i
        end
        def tumor_allelic_depth
          [ tumor_depth - tumor_alt_count, tumor_alt_count ].join ","
        end
      end

      def [] ind
        @muts[ind]
      end

      def each
        @muts.each do |loc,mut|
          yield mut
        end
      end
    end
    class PatchSomaticIndelVcf
      include Pipeline::Task
      requires_file :somaticindel_unpatched_vcf, :somaticindel_verbose
      dumps_file :somaticindel_vcf

      def run
        verbose = SomaticIndelOut.new config.somaticindel_verbose
        unpatched = VCF.read config.somaticindel_unpatched_vcf
        unpatched.preamble_lines.push "##FORMAT=<ID=MM,Number=2,Type=Integer,Description=\"Average number of mismatches across consensus indel and reference-allele supporting reads\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=MQ,Number=2,Type=Integer,Description=\"Average mapping qualitites of consensus indel and reference-allele supporting reads\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=MR,Number=2,Type=Integer,Description=\"Mismatch rate in small 5bp windows around the indel\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=AQ,Number=2,Type=Integer,Description=\"Average base quality computed across all bases in 5bp windows around the indel\">\n"
        unpatched.preamble_lines.push "##FORMAT=<ID=SC,Number=2,Type=Integer,Description=\"Counts of consensus-supporting forward, reverse, reference forward and reverse supporting reads\">\n"
        unpatched.each do |l|
          mut = verbose[ [ l.chrom, l.start.to_i, l.stop.to_i ] ]
          next if !mut
          l.format.push :MM, :MQ, :MR, :AQ, :SC
          l.genotype(config.normal_name).info[:DP] = mut.normal_depth
          l.genotype(config.normal_name).info[:AD] = mut.normal_allelic_depth
          l.genotype(config.normal_name).info[:MM] = mut.n_av_mm.join(",")
          l.genotype(config.normal_name).info[:MQ] = mut.n_av_mapq.join(",")
          l.genotype(config.normal_name).info[:MR] = mut.n_nqs_mm_rate.join(",")
          l.genotype(config.normal_name).info[:AQ] = mut.n_nqs_av_qual.join(",")
          l.genotype(config.normal_name).info[:SC] = mut.n_strand_counts.join(",")

          l.genotype(config.sample_name).info[:DP] = mut.tumor_depth
          l.genotype(config.sample_name).info[:AD] = mut.tumor_allelic_depth
          l.genotype(config.sample_name).info[:MM] = mut.t_av_mm.join(",")
          l.genotype(config.sample_name).info[:MQ] = mut.t_av_mapq.join(",")
          l.genotype(config.sample_name).info[:MR] = mut.t_nqs_mm_rate.join(",")
          l.genotype(config.sample_name).info[:AQ] = mut.t_nqs_av_qual.join(",")
          l.genotype(config.sample_name).info[:SC] = mut.t_strand_counts.join(",")
        end
        unpatched.write config.somaticindel_vcf
      end
    end
  end

  class MutFilter
    include Pipeline::Step
    runs_tasks :filter_muts_pindel
    has_tasks :filter_muts_pindel, :concat_chroms, :filter_muts_annovar, :filter_muts_indelocator, :filter_muts_somatic_indel
    runs_on :tumor_samples, :chroms

    class FilterMuts
      include Pipeline::Task

      EXTRA_HEADERS = [
        :tumor_ref_count, :tumor_alt_count, :tumor_var_freq,
        :normal_ref_count, :normal_alt_count,
        :protein_change, :transcript_change,
        :polyphen2_class, :cosmic_mutations, :segment_logr 
      ]
      ABSOLUTE_HEADERS = [
        :t_ref_count, :t_alt_count, :tumor_var_freq,
        :normal_ref_count, :normal_alt_count, :protein_change,
        :transcript_change, :polyphen2_class, :cosmic_mutations, :segment_logr
      ]

      def mut_to_maf mut
        seg = @segs.find do |seg| 
          seg[:Chromosome] == mut.short_chrom && 
            seg[:Start].to_i < mut.start.to_i &&
            seg[:End].to_i > mut.stop.to_i
        end
        {
          :hugo_symbol => mut.onco.txp_gene, 
          :chromosome => mut.short_chrom,
          :start_position => mut.start,
          :end_position => mut.stop,
          :reference_allele => mut.ref_allele,
          :tumor_seq_allele1 => mut.alt_allele,
          :center => "taylorlab.ucsf.edu",
          :ncbi_build => 37,
          :strand => "+",
          :variant_classification => mut.onco.txp_variant_classification,
          :variant_type => mut.onco.variant_type,
          :dbsnp_rs => (mut.onco.is_snp ? mut.onco.dbSNP_RS : nil), 
          :dbsnp_val_status => mut.onco.dbSNP_Val_Status,
          :tumor_sample_barcode => config.sample_name,
          :matched_norm_sample_barcode => config.normal_name,
          :bam_file => config.sample_bam,
          :protein_change => mut.onco.txp_protein_change,
          :transcript_change => mut.onco.txp_transcript_change,
          :polyphen2_class => mut.onco.pph2_class,
          :cosmic_mutations => mut.onco.Cosmic_overlapping_mutations,
          :segment_logr => seg ? seg[:Segment_Mean].to_f.round(5) : nil
        }
      end
      
      def mutect_to_maf mut
        mut_to_maf(mut).merge({
          :tumor_ref_count => mut.t_ref_count,
          :tumor_alt_count => mut.t_alt_count,
          :t_ref_count => mut.t_ref_count,
          :t_alt_count => mut.t_alt_count,
          :tumor_var_freq => mut.t_var_freq,
          :normal_ref_count => mut.n_ref_count,
          :normal_alt_count => mut.n_alt_count,
        })
      end

      def indel_vcf_to_maf mut
        mut_to_maf(mut).merge({
          :tumor_ref_count => mut.genotype(config.sample_name).ref_count,
          :tumor_alt_count => mut.genotype(config.sample_name).alt_count,
          :t_ref_count => mut.genotype(config.sample_name).ref_count,
          :t_alt_count => mut.genotype(config.sample_name).alt_count,
          :tumor_var_freq => mut.genotype(config.sample_name).alt_freq,
          :normal_ref_count => mut.genotype(config.normal_name).ref_count,
          :normal_alt_count => mut.genotype(config.normal_name).alt_count,
        })
      end

      def create_mafs
        @somatic_maf = Maf.new
        @germline_maf = Maf.new
        @all_muts_maf = Maf.new

        @somatic_maf.headers.concat EXTRA_HEADERS
        @germline_maf.headers.concat EXTRA_HEADERS
        @all_muts_maf.headers.concat ABSOLUTE_HEADERS
        @all_muts_maf.headers.map! do |l|
          l.to_s =~ /_Position/ ? l.to_s.sub(/_Position/,"_position").to_sym : l 
        end
      end

      def write_mafs
        @somatic_maf.sort_by! {|l| -l.tumor_var_freq }
        @germline_maf.sort_by! {|l| -l.tumor_var_freq }
        @somatic_maf.write config.tumor_chrom_maf
        @germline_maf.write config.germline_chrom_maf
        @all_muts_maf.write config.all_muts_chrom_maf
      end

      def load_mutect_snvs chrom
        MuTect.read(config.mutect_snv(chrom), config.mutations_config).each do |l|
          next unless l.keep_somatic? || l.keep_germline?
          log_info "Annotating #{l.contig}:#{l.position}"
          mut = mutect_to_maf l
          unless l.skip_oncotator?
            @somatic_maf.add_line(mut) if l.keep_somatic?
            @germline_maf.add_line(mut) if l.keep_germline? && l.onco.Cosmic_overlapping_mutations
          end
          @all_muts_maf.add_line mut if l.keep_somatic?
        end
      end

      def load_indel_snvs chrom
        VCF.read(indel_vcf(chrom), config.mutations_config).each do |l|
          begin
            log_info "Checking #{l.chrom}:#{l.pos}-#{l.end_pos}"
            next if l.alt.include?("N") || l.ref.include?("N")
            next if l.skip_genotype?([indel_caller, :normal] => config.normal_name) || l.skip_genotype?([indel_caller, :tumor] => config.sample_name)
            next if l.skip_oncotator?
            log_info "Annotating #{l.chrom}:#{l.pos}-#{l.end_pos}"
            mut = indel_vcf_to_maf l
            @somatic_maf.add_line mut
          rescue ArgumentError => e
            log_info e.message
          end
        end
      end

      def run
        create_mafs

        @segs = HashTable.new config.tumor_cnr_seg

        load_mutect_snvs config.chrom
        load_indel_snvs config.chrom

        write_mafs
      end
    end

    class FilterMutsPindel < FilterMuts
      class_init
      requires_files :pindel_vcf, :mutect_snv, :tumor_cnr_seg
      dumps_file :tumor_chrom_maf, :germline_chrom_maf, :all_muts_chrom_maf

      def indel_caller; :pindel; end
      def indel_vcf chrom; config.pindel_vcf chrom; end
    end
    
    class FilterMutsIndelocator < FilterMuts
      class_init
      requires_files :indelocator_vcf, :mutect_snv, :tumor_cnr_seg
      dumps_file :tumor_chrom_maf, :germline_chrom_maf, :all_muts_chrom_maf

      def indel_caller; :indelocator; end
      def indel_vcf chrom; config.indelocator_vcf chrom; end
    end
    class FilterMutsSomaticIndel < FilterMuts
      class_init
      requires_files :somaticindel_vcf, :mutect_snv, :tumor_cnr_seg
      dumps_file :tumor_chrom_maf, :germline_chrom_maf, :all_muts_chrom_maf

      def indel_caller; :somaticindel; end
      def indel_vcf chrom; config.somaticindel_vcf chrom; end
    end

    class ConcatChroms
      include Pipeline::Task
      requires_files :chroms__pindel_vcfs, :chroms__mutect_snvs
      dumps_file :mutect_all_snvs, :pindel_all_vcf

      def run
        mutect = nil
        config.chroms__mutect_snvs.each do |mf|
          m = MuTect.read mf
          if mutect
            mutect.lines.concat m.lines
          else
            mutect = m
          end
        end
        mutect.write config.mutect_all_snvs if mutect
        vcf = nil
        config.chroms__pindel_vcfs.each do |vf|
          v = VCF.read vf
          v.lines.select! do |l|
            (l.ref.size == 1 || l.alt.size == 1) && l.alt !~ /[^ATGC]/  && l.ref !~ /[^ATGC]/
          end
          if vcf
            vcf.lines.concat v.lines
          else
            vcf = v
          end
        end
        vcf.write config.pindel_all_vcf if vcf
      end
    end

    class FilterMutsAnnovar
      include Pipeline::Task
      requires_files :pindel_all_vcf, :mutect_all_snvs
      outs_file :tumor_muts

      def run
        log_info "Filtering mutect and indel output..."
        filter_muts config.mutect_all_snvs, config.pindel_all_vcf, config.tumor_muts or error_exit "Filtering failed"
      end
    end
  end

  class CombineMuts
    include Pipeline::Step
    runs_on :tumor_samples
    runs_tasks :concat_mafs

    class ConcatMafs
      include Pipeline::Task
      requires_files :chroms__tumor_chrom_mafs, :chroms__germline_chrom_mafs, :chroms__all_muts_chrom_mafs
      outs_files :tumor_maf, :germline_maf, :all_muts_maf

      def write_mafs
        @somatic_maf.sort_by! {|l| -l.tumor_var_freq.to_f }
        @germline_maf.sort_by! {|l| -l.tumor_var_freq.to_f }
        @somatic_maf.write config.tumor_maf
        @germline_maf.write config.germline_maf
        @all_muts_maf.write config.all_muts_maf
      end

      def run
        @somatic_maf = Maf.new
        @germline_maf = Maf.new
        @all_muts_maf = Maf.new

        config.sample.chroms.each do |chrom|
          m = Maf.read config.tumor_chrom_maf(chrom)
          @somatic_maf.headers = m.headers
          @somatic_maf.lines.concat m.lines

          m = Maf.read config.germline_chrom_maf(chrom)
          @germline_maf.headers = m.headers
          @germline_maf.lines.concat m.lines

          m = Maf.read config.all_muts_chrom_maf(chrom)
          @all_muts_maf.headers = m.headers
          @all_muts_maf.lines.concat m.lines
        end
        write_mafs
      end
    end
  end
end
