module Pipeline
  module SampleConfig
    extend Pipeline::Config

    def sample(name=nil)
      if name
        samples.find { |s| s[:sample_name] == name }
      else
        samples[sample_index]
      end
    end

    def key_list key
      samples.map { |s| s[key] }.flatten
    end

    def array_key_index_list array, key
      array.map do |s|
        s[key].map_index do |r,i|
          i
        end
      end.flatten
    end

    def array_sample_index_list array, key
      array.map_index do |s,i|
        s[key].map{ i }
      end.flatten
    end

    def sample_index_list key
      array_sample_index_list samples, key
    end

    def key_index_list key
      array_key_index_list samples, key
    end

    def sample_name(name=nil)
      sample(name)[:sample_name]
    end

    def_var :sample_output do |s| "#{output_dir}/#{s || sample_name}" end
  end
end
