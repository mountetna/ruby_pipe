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
  end
end
