module Ribo
  class Babel
    include Pipeline::Step
    runs_tasks :run_babel
    runs_on :babel_tests

    class RunBabel
      include Pipeline::Task
      requires_file :normal_summary
      outs_file :within_babel, :combined_babel, :between_babel

      def run
        r_script :babel, :doBabel, config.normal_summary,
          config.babel_output,
          config.babel.groups.map{|group| "#{group.group_name}=#{group.samples.join(",")}" }.join(":"),
          config.babel_num_reps,
          config.babel_min_rna,
          config.babel_min_rpkm or error_exit "Could not run DESeq"
      end
    end
  end
end
