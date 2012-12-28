require 'extlib'

class Class
  def sister_class(c)
    c = c.to_s.camel_case if c.is_a? Symbol
    (name.split(/::/)[0...-1] + [c]).join("::").to_class
  end

  def daughter_class(c)
    c = c.to_s.camel_case if c.is_a? Symbol
    (name + "::" + c).to_class
  end
end
