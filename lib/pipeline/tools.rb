require 'pipeline/tools/sequence'
require 'pipeline/tools/aligners'
require 'pipeline/tools/mutations'
require 'pipeline/tools/bam'
module Pipeline
  module Tools
    include Pipeline::Tools::Sequence
    include Pipeline::Tools::Aligners
    include Pipeline::Tools::Mutations
    include Pipeline::Tools::Bam

    def cnvkit_fix opts
      args = "#{opts[:target]} #{opts[:background]} #{opts[:reference]} -o #{opts[:output]}"
      cnvkit :fix, args
    end
    def cnvkit_reference opts
      args = format_opts fasta: config.reference_fa, output: opts[:output]
      args += " #{opts[:inputs].join(" ")}"
      cnvkit :reference, args
    end
    def cnvkit_coverage opts
      args = "#{opts[:input]} #{opts[:target]} -o #{opts[:output]}"
      cnvkit :coverage, args
    end
    def cnvkit command, args
      run_cmd "cnvkit.py #{command} #{args}"
    end

    def run_cmd cmd
      log_command cmd
      system "/bin/bash", "-c", cmd
    end
    def read_cmd cmd
      log_command cmd
      output = %x{/bin/bash -c '#{cmd}'}
      $?.success? ? output : nil
    end

    def r_script script, *args
      run_cmd "Rscript #{config.lib_dir}/bin/#{script}.R #{config.lib_dir} #{args.join(" ")}"
    end

    def py_script script, opts
      run_cmd "#{config.lib_dir}/python/#{script}.py #{opts[:args].join(" ")} > #{opts[:out] || '/dev/null'}"
    end

    def java(params)
      opts = {
        :mem => "-Xmx#{params[:mem]}g",
        :tmp => "-Djava.io.tmpdir=#{params[:tmp]}",
        :jar => "-jar #{params[:jar]}",
        :args => params[:args]
      }

      if params[:out]
        run_cmd "#{params[:bin] || config.java} #{ opts.values.join(" ") } > #{params[:out]}"
      else
        run_cmd "#{params[:bin] || config.java} #{ opts.values.join(" ") }"
      end
    end

    def rnaseqc(params)
      java :jar => config.rnaseqc_jar 
    end

    def chimerascan fq1, fq2, output_dir
      run_cmd "python /taylorlab/lib/python/bin/chimerascan_run.py -p #{config.threads} -v --quals #{config.qual_type} #{config.chimerascan_idx} #{fq1} #{fq2} #{output_dir}"
    end

    def cufflinks(opts)
      opts = { :num_threads => config.threads, :quiet => true}.merge(opts)
      bam = opts.delete :bam
      run_cmd "#{config.cufflinks_dir}/cufflinks #{format_opts(opts,true)} #{bam}"
    end

    def cuffmerge(params)
      params = { :gtf => config.reference_gtf, :fa => config.reference_fa, :out => "./merged_asm", :threads => config.threads }.merge(params)
      run_cmd "#{config.cufflinks_dir}/cuffmerge -o #{params[:out]} -g #{params[:gtf]} -s #{params[:fa]} -p #{params[:threads]} #{params[:list]}"
    end

    def cuffcompare(params)
      run_cmd "#{config.cufflinks_dir}/cuffcompare #{params.map{ |o,a| " -#{o} #{a}"}.join(" ") }"
    end

    def cuffdiff(params)
      params = { :gtf => config.reference_gtf, :threads => config.threads }.merge(params)
      run_cmd "#{config.cufflinks_dir}/cuffdiff -p #{params[:threads]} -o #{params[:out]} #{params[:gtf]} #{params[:sample1].join(",")} #{params[:sample2].join(",")}"
    end

    def fastqc(fastq,output)
      run_cmd "fastqc #{fastq} --outdir=#{output}"
    end

    def htmldoc(html,pdf)
      run_cmd "htmldoc --webpage --browserwidth 800px --fontsize 7 -f #{pdf} #{html}" or error_exit "Could not make fastqc pdf"
    end

    def fastx_clipper(params)
      params = { :min_length => 24, :quals => (config.qual_type == "solexa" ? "-Q33" : nil) }.merge(params)
      run_cmd "zcat #{params[:in].join(" ")} | #{config.fastx_dir}/fastx_clipper -a #{params[:adapter]} #{params[:quals]} -l #{params[:min_length]} -n -v -o #{params[:out]}"
    end
    
    def htseq_count(params)
      run_cmd "python -m HTSeq.scripts.count #{params[:input]} #{params[:gtf]} --type=#{params[:type]} --idattr=#{params[:idattr] || :gene_id}  --order=#{params[:order]}  -s no -q > #{params[:out]}"
    end

    def vcf_concat files, outfile
      run_cmd "PERL5LIB=#{config.vcftools_dir}/perl #{config.vcftools_dir}/bin/vcf-concat #{files.join(" ")} > #{outfile}"
    end

    def coverage_bed(bam,intervals,outfile,quality=20)
      run_cmd "samtools view -u -q #{quality} -F 0x400 #{bam} | coverageBed -abam stdin -b #{intervals} -counts > #{outfile}"
    end

    def count_depth bam, chr, pos, pos2=nil
      # get the depth at a particular site
      depth = `samtools mpileup -r #{chr}:#{pos}-#{pos2 || pos} #{bam} 2>/dev/null | awk '{print $4}'`
      return depth.to_i
    end

    def merge_pdfs pdfs, output_pdf
      gs "-o #{output_pdf} -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress #{pdfs.join(" ")}"
    end

    def gs cmd
      run_cmd "#{config.ghostscript_dir}/bin/gs #{cmd}"
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
    def map_opts mmap, opts, sep=" "
      # build an arg string from the map and the method options
      mmap.map do |opt,flag|
        next unless opts[opt]
        opts[opt] = nil if opts[opt] == true

        flag = yield flag if flag

        [ flag, opts[opt] ].compact.join sep
      end.join sep
    end

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
