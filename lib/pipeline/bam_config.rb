module Pipeline
  module BamConfig
    # this just has procs for bam files
    alias_method :load_proc_before_bam, :load_proc
    def load_proc
    end
  end
end
