module Pipeline
  module Tools
    def run_cmd cmd
      log_command cmd
      system cmd
    end

    def r_script script, *args
      run_cmd "#{config.lib_dir}/bin/#{script}.R #{args.join(" ")}"
    end

    def bwa_aln(params)
      params = { :threads => config.threads }.merge(params)
      bwa "aln -t #{params[:threads]} #{config.bwa_idx} #{params[:fq]}", params[:out]
    end

    def bwa_pair(params)
      bwa "sampe #{config.bwa_idx} #{params[:m1]} #{params[:m2]} #{params[:fq1]} #{params[:fq2]}", params[:out]
    end

    def bwa_single(params)
      bwa "samse #{config.bwa_idx} #{params[:map]} #{params[:fq]}", params[:out]
    end

    def bwa(cmd,out)
      run_cmd "#{config.bwa_dir}/bwa #{cmd} > #{out}"
    end

    def picard(jar,params)
      params = { :VALIDATION_STRINGENCY => "SILENT", :TMP_DIR => config.cohort_scratch }.merge(params)
      java :mem => 2, :tmp => config.cohort_scratch, :jar => "#{config.picard_dir}/#{jar.to_s.camel_case}.jar", :out => params.delete(:out), :args => format_opts(params){ |k,v| "#{k}=#{v}" }
    end

    def java(params)
      opts = {
        :mem => "-Xmx#{params[:mem]}g",
        :tmp => "-Djava.io.tmpdir=#{params[:tmp]}",
        :jar => "-jar #{params[:jar]}",
        :args => params[:args]
      }

      if params[:out]
        run_cmd "#{config.java} #{ opts.values.join(" ") } > #{params[:out]}"
      else
        run_cmd "#{config.java} #{ opts.values.join(" ") }"
      end
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
        run_cmd "#{config.samtools_dir}/samtools #{cmd} #{infile} > #{outfile}"
      else
        run_cmd "#{config.samtools_dir}/samtools #{cmd} #{infile}"
      end
    end

    def rnaseqc(params)
      java :jar => config.rnaseqc_jar 
    end

    def gatk(tool,opts)
      opts = { :analysis_type => tool.to_s.camel_case, :reference_sequence => config.reference_fa, :logging_level => "DEBUG", :num_threads => config.threads }.merge(opts)
      java :tmp => config.cohort_scratch, :mem => 4, :jar => "#{config.gatk_dir}/#{config.gatk_jar}", :args => format_opts(opts)
    end

    def mutect(opts)
      opts = { :logging_level => "WARN",
        :analysis_type => "MuTect",
        :baq => "CALCULATE_AS_NECESSARY",
        :reference_sequence => config.reference_fa,
        :dbsnp => config.reference_snp_vcf,
        :num_threads => config.threads,
        :cosmic => config.cosmic_vcf
      }.merge(opts)
      java :mem => 2, :tmp => config.cohort_scratch, :jar => "#{config.mutect_dir}/#{config.mutect_jar}", :args => format_opts(opts)
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

    def tophat(opts)
      opts = { :num_threads => config.threads,
        :output_dir => config.sample_scratch }.merge(opts)
      fq1 = opts.delete :fq1
      fq2 = opts.delete :fq2
      run_cmd "#{config.tophat_dir}/tophat #{format_opts(opts,true)} #{config.bowtie2_idx} #{fq1} #{fq2}"
    end

    def cufflinks(params)
      params = { :threads => config.threads, :gtf => config.reference_gtf }.merge(params)
      gtf = config.find_novel_transcripts ?  nil : "-G #{params[:gtf]}"
      run_cmd "#{config.cufflinks_dir}/cufflinks -q #{gtf} -p #{params[:threads]} -o #{params[:out]} #{params[:bam]}"
    end

    def cuffmerge(params)
      params = { :gtf => config.reference_gtf, :fa => config.reference_fa, :out => "./merged_asm", :threads => config.threads }.merge(params)
      run_cmd "#{config.cufflinks_dir}/cuffmerge -o #{params[:out]} -g #{params[:gtf]} -s #{params[:fa]} -p #{params[:threads]} #{params[:list]}"
    end

    def cuffcompare(params)
      run_cmd "#{config.cufflinks_dir}/cuffcompare #{params.map{ |o,a| " -#{o} #{a}"}.join(" ") }"
    end

    def cuffdiff(params)
      run_cmd "#{config.cufflinks_dir}/cuffdiff -o #{params[:out]} #{params[:gtf]} $BAM1 $BAM2"
    end

    def fastx_clipper(params)
      params = { :min_length => 24, 
        :quals => (config.qual_type == "solexa" ? "-Q33" : nil) }.merge(params)
      if params[:in] =~ /.gz$/
        run_cmd "zcat #{params[:in]} | #{config.fastx_dir}/fastx_clipper -a #{params[:adapter]} #{params[:quals]} -l #{params[:min_length]} -n -v -o #{params[:out]}"
      else
        run_cmd "#{config.fastx_dir}/fastx_clipper -a #{params[:adapter]} #{params[:quals]} -l #{params[:min_length]} -n -v -i #{params[:in]} -o #{params[:out]}"
      end
    end

    def htseq_count(params)
      run_cmd "python -m HTSeq.scripts.count #{params[:input]} #{params[:gtf]} --type=#{params[:type]} -s no -q > #{params[:out]}"
    end

    def filter_muts(snvs,indels,out_file)
      filter_config ="#{config.config_dir}/#{config.genome}_mutationConfig.cfg"
      run_cmd "python #{config.lib_dir}/FilterMutations/Filter.py --keepTmpFiles --tmp #{config.cohort_scratch} #{config.filter_config || filter_config} #{snvs} #{indels} #{out_file}"
    end

    def vcf_concat files, outfile
      run_cmd "PERL5LIB=#{config.vcftools_dir}/perl #{config.vcftools_dir}/bin/vcf-concat #{files.join(" ")} > #{outfile}"
    end

    def hash_table(file,headers=nil)
      lines = File.foreach(file).to_a.map{|s| s.chomp.split(/\t/)}
      headers = lines.shift.map(&:to_sym) if !headers
      lines.map{|l| Hash[headers.zip(l)]}
    end

    def coverage_bed(bam,intervals,outfile)
      run_cmd "#{config.bedtools_dir}/coverageBed -abam #{bam} -b #{intervals} -counts > #{outfile}"
    end

    def count_depth bam, chr, pos, pos2=nil
      # get the depth at a particular site
      depth = `samtools mpileup -r #{chr}:#{pos}-#{pos2 || pos} #{bam} 2>/dev/null | awk '{print $4}'`
      return depth.to_i
    end

    def create_interval_bed
      return if File.exists? config.interval_bed
      File.open(config.interval_bed, "w") do |f|
        File.foreach(config.interval_list) do |g|
          next if g =~ /^\@/
          f.print g
        end
      end
    end

    private
    def format_opts(o,fix_score=nil,&block)
      o.map do |key,value| 
        key = key.to_s.gsub(/_/,"-") if fix_score
        case value
        when Array
          value.map do |v| 
            block ?  block.call(key,v) : "--#{key} #{v}"
          end.join(" ")
        when true
          block ?  block.call(key,nil) : "--#{key}"
        when nil
          nil
        else
          block ?  block.call(key,value) : "--#{key} #{value}"
        end
      end.compact.join(" ")
    end
  end
end


