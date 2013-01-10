class HashTable
  include Enumerable
  def initialize(file,header=nil,comment=nil)
    @lines = File.foreach(file).to_a.map do |s| 
      next if comment && s =~ /^#{comment}/
      s.chomp.split(/\t/)
    end.compact
    @header = header || @lines.shift.map(&:to_sym)
    @lines.map!{|l| Hash[@header.zip(l)]}
  end

  def header
    @header
  end

  def [](ind)
    @lines[ind]
  end

  def print(file)
    File.open(file,"w") do |f|
      f.puts @header.join("\t")
      @lines.each do |l|
        next if l[:_invalid]
        f.puts @header.map{|h| l[h]}.join("\t")
      end
    end
  end

  def each
    @lines.each do |l|
      yield l
    end
  end
end
