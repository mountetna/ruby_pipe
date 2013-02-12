module Pipeline
  module SampleConfig
    extend Pipeline::Config

    def sample(name=nil)
      return samples.find { |s| s[:sample_name] == name } if name
      if job_array
        obj = job_array[job_index]
        while obj
          return obj if obj[:sample_name]
          obj = obj.parent
        end
      end
    end

    def_var :sample_output do |s| "#{output_dir}/#{s || sample_name}" end
  end
end
