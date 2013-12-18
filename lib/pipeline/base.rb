module Pipeline
  module Base
    module ClassMethods
      def class_var *vars
        vars.each do |var|
          # make an accessor
          define_method(var) do
            self.class.send var
          end
        end
      end
    end

    def self.included(base)
      base.extend(ClassMethods)
    end

    def open_error_pid
      ensure_dir config.error_pid

      File.open config.error_pid, "a" do |f|
        yield f
      end
    end

    def ensure_dir(*dirs)
      dirs.each do |f|
        f = File.dirname(f) if !File.directory? f
        FileUtils.mkdir_p f
      end
    end
  end
end
