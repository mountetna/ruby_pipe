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
        set_member key,value
      end.reduce :merge
    end

    def set_member key,value
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
    end

    def add_member key, value
      @obj.update set_member(key,value)
    end

    def owner prop
      if @obj[prop]
        self
      elsif parent && parent.respond_to?(:owner)
        parent.owner(prop) 
      else
        nil
      end
    end

    def property prop
      o = owner(prop)
      o ? o.send(prop) : nil
    end

    def extend_with opts
      key, value = opts.first

      add_member key, value
      @obj[key]
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
      nil
    end
  end
end
