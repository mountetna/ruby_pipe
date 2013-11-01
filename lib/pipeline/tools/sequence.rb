module Pipeline
  module Tools
    module Blat
      def blat_start_server
        run_cmd "gfServer start n0 33333 -canStop -log=#{config.log_dir}/blatServer.log #{config.reference_2bit} &"
      end
      def blat_check_status
        run_cmd "sleep 5"
        run_cmd "gfServer status n0 33333"
      end
      def blat_abort_server
        run_cmd "pkill -9 -f gfServer"
      end

      def blat_stop_server
        run_cmd "gfServer stop n0 33333"
      end
    end
    module Sequence
      include Pipeline::Tools::Blat
    end
  end
end
