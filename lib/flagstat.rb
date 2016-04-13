class Flagstat
  HEADERS= [
    :total, :duplicates, :mapped, :paired_in_sequence, :read1, :read2,
    :properly_paired, :both_mapped, :singletons, :mate_mapped_chr,
    :mate_mapped_chr_highq
  ]

  HEADERS.each do |name|
    define_method(name) do
      @flags[name].first.to_i
    end
    define_method(:"chastity_#{name}") do
      return @flags[name].last.to_i
    end
  end
  
  def initialize file
    @flags = Hash[
      HEADERS.zip(
        File.foreach(file).each_with_index.map do |l|
          l.scan(/([0-9]+) \+ ([0-9]+)/).flatten
        end
      )
    ]
  end

  def each
    @flags.each do |flag,l|
      yield f,l
    end
  end

  def clean_hash
    HEADERS.map do |name|
      { name => send(name) }
    end.reduce(:merge)
  end

  def chastity_hash
    HEADERS.map do |name|
      { name => send(:"chastity_#{name}") }
    end.reduce(:merge)
  end

  def method_missing(method, *args, &block)
    return @flags[method].first.to_i if @flags[method]
    method.to_s.match(/^chastity_(.*)/) do |m|
      return @flags[m[1].to_sym].last.to_i if @flags[m[1].to_sym]
    end
    super
  end
end
