require 'pipeline'
require 'exome/align'
require 'exome/merge'
require 'exome/recal'
require 'exome/hybrid_qc'
require 'exome/mut_det'
require 'exome/config'
require 'exome/copy_number'

module Exome
  class PairedAlign 
    include Pipeline::Script
    runs_steps :align, :merge, :recal, :hybrid_qc, :mut_det, :copy_number

    def exclude_task? task
      case task
      when Exome::CopyNumber::CopySeg, Exome::CopyNumber::MutAddSeg
        return true if config.unit != "exome"
      when Exome::CopyNumber::MutAddCnv
        return true if config.unit == "exome"
      end
      return nil
    end

    class ConfigGenerator
      include Pipeline::ConfigGenerator

      def input(args)
        sample = get_sample args.first
        sample[:inputs] = []
        getlines "your fastq files" do |l|
          sample[:inputs] += l.strip.split
        end
        sample[:inputs] = pair_fastqs sample[:inputs]
      end
      usage "input <sample name>", "Input fastqs for the given sample"

      def intervals(args)
        config[:interval_list] = args.first
      end
      usage "intervals <interval list>", "Set the interval list file for the given sample."

      def normal(args)
        sample = find_sample args[0]
        normal = find_sample args[1]
        if !sample || !normal
          puts "No sample of that name.".red
          return
        end
        if sample == normal
          puts "Sample and normal are the same".red
          return
        end

        sample[:normal_name] = normal[:sample_name]
      end
      usage "normal <sample name> <normal sample name>", "Set the normal name for the given sample (the first sample is assumed to be the normal)"
    end

    module Config
      include Pipeline::SampleConfig
      # this is config stuff that is particular to this sample

      def sample_index
        # just tell the index of the current sample
        case step
        when :align
          sample_index_list(:input_fastq1_list)[job_index]
        when :merge, :hybrid_qc
          job_index
        when :mut_det, :copy_number
          job_index+1
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
        when :mut_det, :copy_number
          # this runs once on every tumor sample
          samples.size - 1
        end
      end

      def scratch
        case step
        when :align, :merge, :hybrid_qc, :mut_det, :copy_number
          sample_scratch
        else
          job_scratch
        end
      end
    end
  end
end
