module Pipeline
  module Tools
    module IndelDetection
      def indelocator(opts)
        opts = { :logging_level => config.logging_level,
          :baq => "CALCULATE_AS_NECESSARY",
          :reference_sequence => config.reference_fa,
          :analysis_type => "IndelGenotyperV2"
        }.merge(opts)
        java :mem => 2, :tmp => config.sample_tmp, :jar => "#{config.indelocator_dir}/#{config.indelocator_jar}", :args => format_opts(opts)
      end

      def somatic_indel_detector(opts)
        opts = { :logging_level => config.logging_level,
          :baq => "CALCULATE_AS_NECESSARY",
          :reference_sequence => config.reference_fa,
          :analysis_type => :SomaticIndelDetector
        }.merge(opts)
        java :mem => 2, :tmp => config.sample_tmp, :jar => "#{config.somaticindel_dir}/#{config.somaticindel_jar}", :args => format_opts(opts)
      end

      def pindel(opts)
        opts = { :number_of_threads => config.threads, :fasta => config.reference_fa }.merge(opts)

        tempfile = opts.delete :tempfile
        File.open(tempfile,"w") do |f| 
          opts[:bams].each do |b|
            f.puts "#{b[:bam]} #{config.frag_size} #{b[:name]}"
          end
        end
        opts[:"config-file"] = tempfile
        opts.delete :bams

        run_cmd "#{config.pindel_dir}/pindel #{format_opts(opts)}"
      ensure
        File.unlink(tempfile) if tempfile && File.exists?(tempfile)
      end

      def pindel_to_vcf(opts)
        opts = { :reference => config.reference_fa, :reference_name => config.reference_name, :reference_date => config.reference_date }.merge(opts)
        run_cmd "#{config.pindel_dir}/pindel2vcf #{format_opts(opts)}"
      end

      def strelka(params)
        params = { :reference => config.reference_fa }.merge(params)

        #require 'fasta'
        #f = Fasta.new params[:reference]
        run_cmd "#{config.strelka_dir}/libexec/strelka2 -bam-file #{params[:normal]} --tumor-bam-file #{params[:tumor]} -clobber --somatic-indel-file #{params[:output]} -bam-seq-name #{params[:chrom]} -samtools-reference #{params[:reference]} -genome-size #{3095693984} --tumor-realigned-read-file #{params[:realigned_bam]}"
      end
    end
    module SNVDetection
      def mutect(opts)
        opts = { 
          :logging_level => config.logging_level,
          :analysis_type => "MuTect",
          :baq => "CALCULATE_AS_NECESSARY",
          :reference_sequence => config.reference_fa,
          :dbsnp => config.reference_snp_vcf,
          :num_threads => config.threads,
          :cosmic => config.cosmic_vcf,
        }.merge(opts)

        opts = {
          :max_alt_alleles_in_normal_count => 100000,
          :max_alt_alleles_in_normal_qscore_sum => 100000,
          :max_alt_allele_in_normal_fraction => 1
        }.merge(opts) if opts.delete :no_normal_filter

        java :bin => config.java6, :mem => 2, :tmp => config.sample_tmp, :jar => "#{config.mutect_dir}/#{config.mutect_jar}", :args => format_opts(opts)
      end

      def freebayes(opts)
        opts = { 
          :fasta_reference => config.reference_fa,
        }.merge(opts)
        bam_file = opts.delete :bam
        output_file = opts.delete :output
        run_cmd "#{config.freebayes_dir}/bin/freebayes #{format_opts(opts,true)} #{bam_file} > #{output_file}"
      end
    end
    module Annotation
      def filter_muts(snvs,indels,out_file)
        filter_config ="#{config.config_dir}/#{config.genome}_mutationConfig.cfg"
        run_cmd "python #{config.lib_dir}/FilterMutations/Filter.py --keepTmpFiles --tmp #{config.sample_tmp} #{config.filter_config || filter_config} #{snvs} #{indels} #{out_file}"
      end

      def snpeff infile, outfile
        args = [
          "eff",
          "-canon",
          "-noLog",
          "-t",
          "-v GRCh37.75",
          infile ].join(" ")
        java :mem => 2, :tmp => config.sample_tmp, :jar => "#{config.snpeff_dir}/#{config.snpeff_jar}", :args => args, :out => outfile
      end
    end
    module Mutations
      include Pipeline::Tools::IndelDetection
      include Pipeline::Tools::SNVDetection
      include Pipeline::Tools::Annotation
    end
  end
end
