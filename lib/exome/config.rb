module Exome
  class Config
    include Pipeline::Config

    def input_fastq1_list
      CONFIG['INPUT_FASTQ1_LIST'].split(/,/).map{|c| c.split(/:/)}
    end

    def input_fastq2_list
      CONFIG['INPUT_FASTQ2_LIST'].split(/,/).map{|c| c.split(/:/)}
    end

    def input_fastq1
      # just take 'em one at a time.
    end
  end
end
