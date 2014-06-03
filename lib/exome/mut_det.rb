#!/usr/bin/env ruby
require 'hash_table'
require 'fileutils'
require 'mutect'
require 'vcf'
require 'maf'

module Exome
  class MutDet
    include Pipeline::Step
    runs_tasks :mutect, :pindel, :pindel_vcf, :patch_pindel_vcf
    has_tasks :mutect, :pindel, :pindel_vcf, :patch_pindel_vcf, :strelka, :somatic_indel_detector, :patch_somatic_indel_vcf
    resources :threads => 1
    runs_on :tumor_samples, :chroms

    def vacuum
      config.sample_names.each do |s|
        FileUtils.rm Dir.glob("#{config.scratch_dir}/#{s}/AnnoVar*")
      end
    end

    class Mutect
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_files :mutect_snvs, :mutect_coverage

      def run
	      log_info "Running muTect for tumor #{config.sample_name}, normal #{config.normal_name}"
        mutect "input_file:normal" => config.normal_bam, "input_file:tumor" => config.tumor_bam,
          :intervals => config.chrom.chrom_name,
          :no_normal_filter => config.disable_mutect_normal_filter,
          :out => config.mutect_snvs_tmp, :coverage_file => config.mutect_coverage or error_exit "muTect failed"

        # kludge to make sure mutect completes before ensuring this step
        FileUtils.mv config.mutect_snvs_tmp, config.mutect_snvs
      end
    end

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

    class PindelVcf
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
        unpatched.each do |l|
          # skip invalid ones
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
         :n_obs_counts => [ :alt, :any, :total],
         :n_av_mm => [:alt,:ref],
         :n_av_mapq => [:alt,:ref],
         :n_nqs_mm_rate => [:alt,:ref],
         :n_nqs_av_qual => [:alt,:ref],
         :n_strand_counts => [:alt_forward,:alt_reverse,:ref_forward,:ref_reverse],
         :t_obs_counts => [:alt,:any,:total],
         :t_av_mm => [:alt,:ref],
         :t_av_mapq => [:alt,:ref],
         :t_nqs_mm_rate => [:alt,:ref],
         :t_nqs_av_qual => [:alt,:ref],
         :t_strand_counts => [:alt_forward,:alt_reverse,:ref_forward,:ref_reverse],
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

        def normal_nonref_count
          @normal_nonref_count ||= n_obs_counts[:any].to_i
        end

        def normal_allelic_depth
          [ normal_depth - normal_nonref_count, normal_alt_count ].join ","
        end

        def tumor_depth
          @tumor_depth ||= t_obs_counts[:total].to_i
        end

        def tumor_alt_count
          @tumor_alt_count ||= t_obs_counts[:alt].to_i
        end

        def tumor_nonref_count
          @tumor_nonref_count ||= t_obs_counts[:any].to_i
        end

        def tumor_allelic_depth
          [ tumor_depth - tumor_nonref_count, tumor_alt_count ].join ","
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
        unpatched.each do |l|
          mut = verbose[ [ l.chrom, l.start.to_i, l.stop.to_i ] ]
          next if !mut
          l.genotype(config.normal_name).info[:DP] = mut.normal_depth
          l.genotype(config.normal_name).info[:AD] = mut.normal_allelic_depth

          l.genotype(config.sample_name).info[:DP] = mut.tumor_depth
          l.genotype(config.sample_name).info[:AD] = mut.tumor_allelic_depth
        end
        unpatched.write config.somaticindel_vcf
      end
    end

    class Strelka
      include Pipeline::Task
      requires_files :normal_bam, :tumor_bam
      dumps_file :strelka_vcf, :strelka_bam
      
      def run
        strelka :tumor => config.tumor_bam, :normal => config.normal_bam, :output => config.strelka_vcf, :chrom => config.chrom.chrom_name, :realigned_bam => config.strelka_bam or error_exit "Strelka failed"
      end
    end
  end

  class MutFilter
    include Pipeline::Step
    runs_tasks :filter_muts_pindel
    has_tasks :filter_muts_pindel, :concat_chroms, :filter_muts_annovar, :filter_muts_somatic_indel
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
        MuTect.read(config.mutect_snvs(chrom), config.mutations_config).each do |l|
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
      requires_files :pindel_vcf, :mutect_snvs, :tumor_cnr_seg
      dumps_file :tumor_chrom_maf, :germline_chrom_maf, :all_muts_chrom_maf

      def indel_caller; :pindel; end
      def indel_vcf chrom; config.pindel_vcf chrom; end
    end

    class FilterMutsSomaticIndel < FilterMuts
      class_init
      requires_files :somaticindel_vcf, :mutect_snvs, :tumor_cnr_seg
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
