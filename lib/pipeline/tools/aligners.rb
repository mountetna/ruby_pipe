module Pipeline
  module Tools
    module BWA
      def bwa_aln(params)
        params = { :threads => config.threads }.merge(params)
        if params[:bam]
          bwa "aln -t #{params[:threads]} #{config.bwa_idx} -b#{params[:pair]} #{params[:bam]}", params[:out]
        else
          bwa "aln -t #{params[:threads]} #{config.bwa_idx} #{params[:fq]}", params[:out]
        end
      end

      def bwa_pair(params)
        bwa "sampe #{config.bwa_idx} #{params[:m1]} #{params[:m2]} #{params[:fq1]} #{params[:fq2]}", params[:out]
      end


      def bwa_mem(params)
        params = { threads: config.threads, index: config.bwa_idx, mark_secondary: true }.merge(params)

        if params[:fq1].is_a?(Array)
          params[:fq1] = "<(zcat #{params[:fq1].join(" ")})"
        end
        if params[:fq2].is_a?(Array)
          params[:fq2] = "<(zcat #{params[:fq2].join(" ")})"
        end
        opt_map = { threads: :t,
                    mark_secondary: :M,
                    mismatch_penalty: :B,
                    min_score: :T,
                    index: nil,
                    fq1: nil,
                    fq2: nil }
        out = params.delete :out
        args = map_opts(opt_map, params) { |flag| "-#{flag}" }
        bwa "mem #{args}", out
      end

      def bwa_single(params)
        bwa "samse #{config.bwa_idx} #{params[:map]} #{params[:fq]}", params[:out]
      end

      def bwa(cmd,out)
        run_cmd "#{config.bwa_dir}/bwa #{cmd} > #{out}"
      end
    end
    module Tophat
      def tophat(opts)
        opts = { :num_threads => config.threads }.merge(opts)
        fq1 = opts.delete :fq1
        fq2 = opts.delete :fq2
        run_cmd "#{config.tophat_dir}/tophat #{format_opts(opts,true)} #{config.bowtie2_idx} #{fq1} #{fq2}"
      end
    end
    module Rsem
      def rsem_format *list
        list.map do |s| 
          file = File.realdirpath(s)
          file =~ /.gz$/ ? "<(zcat #{file})" : file
        end.join(",")
      end

      def rsem type, opts
        #opts[ "#{config.qual_type}_quals".to_sym ] = true if !opts[:bam] && !opts[:sam]
        args = opts.delete(:args)
        output = opts.delete :output
        ensure_dir opts[:temporary_folder]

        Dir.chdir(output) do
          run_cmd "#{config.rsem_dir}/rsem-#{type.to_s.gsub(/_/,"-")} #{format_opts(opts,true)} #{args.values.join(" ")}"
        end
      end
    end
    module Kallisto
      def kallisto_quant opts
        fastqs = opts.delete :input
        bam = opts.delete :bam
        opts = opts.merge(
          pseudobam: true
        )
        run_cmd "#{config.kallisto_dir}/kallisto quant #{format_opts(opts,true)} #{fastqs.join(" ")} | samtools view -Sb - > #{bam}"
      end
    end
    module Aligners
      include Pipeline::Tools::BWA
      include Pipeline::Tools::Tophat
      include Pipeline::Tools::Rsem
      include Pipeline::Tools::Kallisto
    end
  end
end
