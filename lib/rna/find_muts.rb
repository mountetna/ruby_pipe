require 'pipeline'
require 'rna/univ_geno'
require 'rna/config'

module Rna
  class FindMuts 
    include Pipeline::Script
    runs_steps :univ_geno

    module Config
      # this is config stuff that is particular to this sample
      def sample_name
        case step
        when :univ_geno
          job_name
        end
      end

      def splits
        case step
        when :univ_geno
          # this runs once
          1
        end
      end

      def scratch
        case step
        when :univ_geno
          job_scratch
        end
      end
    end
  end
end
