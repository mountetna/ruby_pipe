module Pipeline
  module Tools
    def bwa_aln(params)
      params = { :threads => config.resources[:threads] }.merge(params)
      bwa "aln -t #{params[:threads]} #{config.bwa_idx} #{params[:fq]}", params[:out]
    end

    def bwa_pair(params)
      bwa "sampe #{config.bwa_idx} #{params[:m1]} #{params[:m2]} #{params[:fq1]} #{params[:fq2]}", params[:out]
    end

    def bwa(cmd,out)
      system "#{config.bwa_dir}/bwa #{cmd} > #{out}"
    end

    def picard(jar,params)
      params = { :VALIDATION_STRINGENCY => "SILENT", :TMP_DIR => config.cohort_scratch }.merge(params)
      java :mem => 2, :tmp => config.cohort_scratch, :jar => "#{config.picard_dir}/#{jar.to_s.camel_case}.jar", :args => format_opts(params){ |k,v| "#{k}=#{v}" }
    end

    def java(params)
      opts = {
        :mem => "-Xmx#{params[:mem]}g",
        :tmp => "-Djava.io.tmpdir=#{params[:tmp]}",
        :jar => "-jar #{params[:jar]}",
        :args => params[:args]
      }
      params.keys.each do |p|
        params[p] = opts[p]
      end

      system "#{config.java} #{ params.values.join(" ") }"
    end

    def sam_to_bam(infile,outfile)
      samtools "view -bS", infile, outfile
    end

    def sam_sample_name(bam)
      IO.popen("#{config.samtools_dir}/samtools view -H #{bam}").readlines.map{ |l| l.match(/^@RG.*SM:(.*?)\s/) { |m| m[1] } }.compact.first
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
        system "#{config.samtools_dir}/samtools #{cmd} #{infile} > #{outfile}"
      else
        system "#{config.samtools_dir}/samtools #{cmd} #{infile}"
      end
    end

    def rnaseqc(params)
      java :jar => config.rnaseqc_jar 
    end

    def gatk(tool,opts)
      opts = { :analysis_type => tool.to_s.camel_case, :reference_sequence => config.hg19_fa, :logging_level => "DEBUG" }.merge(opts)
      java :tmp => config.cohort_scratch, :mem => 4, :jar => "#{config.gatk_dir}/#{config.gatk_jar}", :args => format_opts(opts){ |k,v| "--#{k} #{v}" }
    end

    def mutect(opts)
      opts = { :logging_level => "WARN", :analysis_type => "MuTect", :baq => "CALCULATE_AS_NECESSARY", :reference_sequence => config.hg19_fa }.merge(opts)
      java :mem => 2, :tmp => config.cohort_scratch, :jar => "#{config.mutect_dir}/#{config.mutect_jar}", :args => format_opts(opts){ |k,v| "--#{k} #{v}" }
    end

    def tophat(params)
      params = { :threads => config.resources[:threads], :scratch => config.cohort_scratch, :gtf => config.hg19_ucsc_gtf }.merge(params)
      system "#{config.tophat_dir}/tophat -G #{params[:gtf]} -o #{params[:scratch]} -r #{params[:frag_size]} -p #{params[:threads]} #{config.bowtie2_idx} #{params[:fq1]} #{params[:fq2]}"
    end

    def cufflinks(params)
      params = { :threads => config.resources[:threads], :gtf => config.hg19_ucsc_gtf }.merge(params)
      system "#{config.cufflinks_dir}/cufflinks -q -p #{params[:threads]} -o #{params[:out]} #{params[:bam]}"
    end

    def cuffmerge(params)
      params = { :gtf => config.hg19_ucsc_gtf, :fa => config.hg19_fa, :out => "./merged_asm", :threads => config.resources[:threads] }.merge(params)
      system "#{config.cufflinks_dir}/cuffmerge -o #{params[:out]} -g #{params[:gtf]} -s #{params[:fa]} -p #{params[:threads]} #{params[:list]}"
    end

    def cuffcompare(params)
      system "#{config.cufflinks_dir}/cuffcompare #{params.map{ |o,a| " -#{o} #{a}"}.join(" ") }"
    end

    def cuffdiff(params)
      system "#{config.cufflinks_dir}/cuffdiff -o #{params[:out]} #{params[:gtf]} $BAM1 $BAM2"
    end

    def fastx_clipper(args)
      system "#{config.fastx_dir}/fastx_clipper #{args}"
    end

    def htseq_count(args)
      system "python -m HTSeq.scripts.count #{ args }"
    end

    def filter_muts(snvs,indels,out_file)
      filter_config ="#{config.lib_dir}/FilterMutations/mutationConfig.cfg"
      system "python #{config.lib_dir}/FilterMutations/Filter.py --keepTmpFiles --tmp #{config.cohort_scratch} #{config.filter_config || filter_config} #{snvs} #{indels} #{out_file}"
    end

    def hash_table(file,headers=nil)
      lines = File.foreach(file).to_a.map{|s| s.chomp.split(/\t/)}
      headers = lines.shift.map(&:to_sym) if !headers
      lines.map{|l| Hash[headers.zip(l)]}
    end

    def coverage_bed(bam,intervals,outfile)
      system "#{config.bedtools_dir}/coverageBed -abam #{bam} -b #{intervals} -counts > #{outfile}"
    end

    private
    def format_opts(o,&block)
      o.map do |key,value| 
        value.is_a?(Array) ? value.map do |v| 
          yield key,v
        end.join(" ") : yield(key,value)
      end.join(" ")
    end
  end
end


