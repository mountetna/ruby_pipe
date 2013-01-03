require 'pipeline'
require 'exome/align'
require 'exome/merge'
require 'exome/recal'
require 'exome/hybrid_qc'
require 'exome/mut_det'
require 'exome/config'
#require 'exome/copy_number'

module Exome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :align, :merge, :recal, :hybrid_qc, :mut_det

    module Config
      # this is config stuff that is particular to this sample
      def sample_name
        case step
        when :align
          samples[samples.to_enum.with_index.map{ |s,i| s[:input_fastq1_list].map{ i } }.flatten[job_index]][:sample_name]
        when :merge, :hybrid_qc
          samples[job_index][:sample_name]
        when :mut_det
          samples[job_index+1][:sample_name]
        end
      end

      def splits
        case step
        when :align
          # this runs once on every input fastq
          input_fastq1_list.size
        when :merge, :hybrid_qc
          # this runs once on every sample
          sample_names.size
        when :recal
          # this runs once
          1
        when :mut_det
          # this runs once on every tumor sample
          sample_names.size - 1
        end
      end

      def scratch
        case step
        when :align, :merge, :hybrid_qc, :mut_det
          sample_scratch
        else
          job_scratch
        end
      end
    end
  end
end
