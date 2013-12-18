#!/usr/bin/env ruby

module Exome
  class StartBlat
    include Pipeline::Step
    runs_tasks :start_server
    runs_on :cohort
    resources :node => 0

    def run
      blat_start_server

      seconds = 0
      until blat_check_status
        seconds += 5
        if seconds > 120
          blat_abort_server
          error_exit "Couldn't talk to blat server"
        end
      end
    end
  end
  class DetectRearrangements
    include Pipeline::Step
    runs_tasks :start_blat, :extract_sclips, :run_crest, :stop_blat
    runs_on :tumor_samples
    audit_report :sample_name, :normal_name
  end
end
