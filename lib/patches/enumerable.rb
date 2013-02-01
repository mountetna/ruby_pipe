module Enumerable
  def map_index &block
    to_enum.with_index.map do |o,i|
      yield o, i
    end
  end
end
