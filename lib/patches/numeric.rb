class Numeric
  def to_human
    # make this human-readable
    return self.to_s if self < 1000
    [ :T, :G, :M, :K ].each_with_index do |s,i|
      n = self.to_f / 10**((4-i)*3)
      if n > 1
        return "#{n.round(1)}#{s}"
      end
    end
    return self.to_s
  end
end
