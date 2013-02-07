module Pipeline
  class SampleObject
    attr_reader :parent, :index
    def initialize(o,p,i=nil)
      @parent = p
      @obj = replace_members o
      @index = i
    end

    def replace_members o
      o.map do |key,value|
        {
          key => case value
          when Hash
            SampleObject.new(value,self)
          when Array
            value.map_index do |h,i|
              if h.is_a? Hash
                SampleObject.new(h,self,i)
              else
                h
              end
            end
          else
            value
          end
        }
      end.reduce :merge
    end

    def [](key)
      @obj[key]
    end

    def []=(key,value)
      @obj[key] = value
    end

    def method_missing(meth,*args,&block)
      if @obj[meth.to_sym]
        return @obj[meth.to_sym]
      end
      super
    end
  end
end
