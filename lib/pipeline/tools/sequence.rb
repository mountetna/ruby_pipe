module Pipeline
  module Tools
    module Blat
      def blat_start_server 
        run_cmd "#{config.blat_dir}/gfServer start #{config.blat_node} #{config.blat_port} -canStop -log=#{config.log_dir}/blatServer.log #{config.reference_2bit}"
      end
      def blat_check_status
        run_cmd "sleep 5"
        run_cmd "#{config.blat_dir}/gfServer status #{config.blat_node} #{config.blat_port}"
      end
      def blat_abort_server
        run_cmd "pkill -9 -f gfServer"
      end

      def blat_stop_server
        run_cmd "#{config.blat_dir}/gfServer stop #{config.blat_node} #{config.blat_port}"
      end
    end
    module Rearrangement
      def extract_sclips opts
        opts = { :dir => config.sample_scratch, :sample => config.sample_name }.merge(opts)
        run_cmd "#{config.crest_dir}/extractSClip.pl -i #{opts[:bam]} --ref_genome #{config.reference_fa} -p #{opts[:sample]} -o #{File.realdirpath opts[:dir]} -r #{opts[:chrom]}" 
      end

      def crest opts
        opts = { :dir => File.realdirpath(config.sample_scratch) }.merge(opts)
        run_cmd "#{config.crest_dir}/CREST.pl -f #{opts[:sclip]} -d #{opts[:tumor_bam]} -g #{opts[:normal_bam]} --ref_genome #{config.reference_fa} -t #{config.reference_2bit} -r #{opts[:chrom]} --blatserver #{config.blat_node} --blatport #{config.blat_port} --sensitive --out_dir #{opts[:dir]} --prefix #{opts[:prefix]}"
      end
    end
    module Sequence
      include Pipeline::Tools::Blat
      include Pipeline::Tools::Rearrangement
    end
  end
end
