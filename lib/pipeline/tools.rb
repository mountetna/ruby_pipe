module Pipeline
  module Tools
    def bwa_aln(params)
      params = { :threads => 12 }.merge(params)
      bwa "aln -t #{params[:threads]} #{config.bwa_index} #{params[:fq]}", params[:out]
    end

    def bwa_sampe(params)
      bwa "sampe #{config.bwa_index} #{params[:m1]} #{params[:m2]} #{params[:fq1]} #{params[:fq2]}", params[:out]
    end

    def bwa(cmd,out)
      system "#{config.bwa_dir}/bwa #{cmd} > #{out}"
    end

    def picard(jar,params)
      params = { :VALIDATION_STRINGENCY => "SILENT", :TMP_DIR => config.scratch }.merge(params)
      java :mem => 2, :tmp => config.scratch, :jar => "#{config.picard_dir}/#{jar.to_s.camel_case}.jar", :args => (params.map{ |key,value| value.is_a?(Array) ? value.map{|v| "#{key}=#{v}"}.join(" ") : "#{key}=#{value}" }.join(" "))
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

    def sam_flag(infile,outfile)
      samtools "flagstat", infile, outfile
    end

    def sam_index(infile)
      samtools "index", infile, nil
    end

    def samtools(cmd,infile,outfile)
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
      opts = { "analysis_type" => tool.to_s.camel_case, "reference_sequence" => config.hg19_fa, "logging_level" => "DEBUG" }.merge(opts)
      java :tmp => config.scratch, :mem => 4, :jar => "#{config.gatk_dir}/#{config.gatk_jar}", :args => opts.map { |flag,opt| "--#{flag} #{opt}" }.join(" ")
    end

    def mutect(opts)
      opts = { "logging-level" => "WARN", "analysis_type" => "MuTect", "baq" => "CALCULATE_AS_NECESSARY", "reference-sequence" => config.hg19_fa }.merge(opts)
      java :mem => 2, :tmp => config.scratch, :jar => "#{config.mutect_dir}/#{config.mutect_jar}", :args => opts.map { |flag,opt| "--#{flag} #{opt}" }.join(" ")
    end

    def tophat(args)
      system "#{config.tophat_dir}/tophat #{args}"
    end

    def fastx_clipper(args)
      system "#{config.fastx_dir}/fastx_clipper #{args}"
    end

    def htseq_count(args)
      system "python -m HTSeq.scripts.count #{ args }"
    end

    def filter_muts(snvs,indels,out_file)
      filter_config ="#{config.lib_dir}/FilterMutations/mutationConfig.cfg"
      system "python #{config.lib_dir}/FilterMutations/Filter.py --keepTmpFiles --tmp #{config.scratch} #{config.filter_config || filter_config} #{snvs} #{indels} #{out_file}"
    end

    def replace_dict(args)
    end

    def coverage_bed(bam,intervals,outfile)
      system "#{config.bedtools_dir}/coverageBed -abam #{bam} -b #{intervals} -counts | cut -f 1-3,5 > #{outfile}"
    end
  end
end
