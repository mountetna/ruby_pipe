require 'pipeline'
require 'exome/config'
require 'exome/copy_number'

module Exome
  class Cnv 
    include Pipeline::Script
    runs_steps :copy_number
    #exclude_task :copy_number, :update_mutations

    def exclude_task? task
      case task
      when Exome::CopyNumber::CopySeg
        return true if config.unit != "exome"
      when Exome::CopyNumber::MutAddSeg, Exome::CopyNumber::MutAddCnv
        return true
      end
      return nil
    end

    module Config
      # this is config stuff that is particular to this sample
      def sample_index
        case step
        when :copy_number
          # list 
          job_index+1
        end
      end

      def splits
        case step
        when :copy_number
          # this runs once on every tumor sample that has a normal bam
          sample_names.size - 1
        end
      end

      def scratch
        case step
        when :copy_number
          sample_scratch
        else
          job_scratch
        end
      end
    end
  end
end
