require 'extlib'

class Class
  def class_chain
    name.split(/::/)
  end

  def parent_chain
    class_chain[0...-1]
  end

  def sister_class(c)
    c = c.to_s.camel_case if c.is_a? Symbol
    join_chain parent_chain + [c]
  end

  def daughter_class(c)
    c = c.to_s.camel_case if c.is_a? Symbol
    join_chain class_chain + [c]
  end

  def parent_class
    join_chain parent_chain
  end

  def class_symbol
    class_chain.map(&:snake_case).join('_').to_sym
  end

  private
  def join_chain chain
    chain.inject(Object) { |m,c| m.const_get(c) }
  end
end
