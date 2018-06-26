module Pipeline
  module Tools
    module Gatk
      def gatk(tool,opts)
        opts = { :analysis_type => tool.to_s.camel_case, :reference_sequence => config.reference_fa, :logging_level => config.logging_level, :num_threads => config.threads }.merge(opts)
        java :tmp => config.sample_tmp, :mem => 4, :jar => "#{config.gatk_dir}/#{config.gatk_jar}", :args => format_opts(opts)
      end
    end
    module Picard
      def picard(jar,params)
        params = { :VALIDATION_STRINGENCY => "SILENT", :TMP_DIR => config.sample_tmp }.merge(params)
        java :mem => 4, :tmp => config.sample_tmp, :jar => "#{config.picard_dir}/#{jar.to_s.camel_case}.jar", :out => params.delete(:out), :args => format_opts(params){ |k,v| "#{k}=#{v}" }
      end
    end
    module Samtools
      def sam_to_bam(infile,outfile)
        samtools "view -bS", infile, outfile
      end

      def sam_merge(outfile, *infiles)
        samtools :merge, "#{outfile} #{infiles.join(" ")}"
      end

      def sam_sample_name(bam)
        IO.popen("#{config.samtools_dir}/samtools view -H #{bam}").readlines.map do |l| 
          m = l.match(/^@RG.*SM:(?<sample_name>.*?)\s/)
          m ? m[:sample_name] : nil
        end.compact.first
      end

      def sam_flags(infile,outfile)
        samtools "flagstat", infile, outfile
      end

      def sam_sort(infile,outpref)
        samtools "sort #{infile}", outpref
      end

      def sam_reads(bam,chr,start,stop)
        `#{config.samtools_dir}/samtools view #{bam} #{chr}:#{start}-#{stop}`.split(/\n/).map do |l|
          # these 
          Hash[
            [ :qname, :flag, :contig, :pos, :mapq, :cigar, :mate_contig, :mate_pos, :insert_size, :seq, :quality, :opt ].zip l.split(/\t/,12)
          ]
        end
      end

      def sam_index(infile)
        samtools "index", infile
      end

      def samtools(cmd,infile,outfile=nil)
        if outfile
          run_cmd "#{config.samtools_dir}/samtools #{cmd} #{infile} > #{outfile}"
        else
          run_cmd "#{config.samtools_dir}/samtools #{cmd} #{infile}"
        end
      end
    end
    module Bedtools
      def bedtools_bamtofastq bam, fq1, fq2
        run_cmd "#{config.bedtools_dir}/bedtools bamtofastq -i #{bam} -fq >(gzip -c > #{fq1}) -fq2 >(gzip -c > #{fq2})"
      end
    end
    module Bam
      include Pipeline::Tools::Gatk
      include Pipeline::Tools::Picard
      include Pipeline::Tools::Samtools
      include Pipeline::Tools::Bedtools
    end
  end
end
