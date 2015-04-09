class SomaticIndelOut
  include Enumerable
  FIELDS = %r{
    (?<field_name> [A-Z_]+ ){0}
    (?<params> [A-Z\/]+ ){0}
    (?<fields> [\d\.\/-]+){0}
    \A\g<field_name>(\[\g<params>\])?:?\g<fields>?\Z
  }x
  FIELD_KEYS = {
     :n_obs_counts => [ :alt, :any, :total],
     :n_av_mm => [:alt,:ref],
     :n_av_mapq => [:alt,:ref],
     :n_nqs_mm_rate => [:alt,:ref],
     :n_nqs_av_qual => [:alt,:ref],
     :n_strand_counts => [:alt_forward,:alt_reverse,:ref_forward,:ref_reverse],
     :t_obs_counts => [:alt,:any,:total],
     :t_av_mm => [:alt,:ref],
     :t_av_mapq => [:alt,:ref],
     :t_nqs_mm_rate => [:alt,:ref],
     :t_nqs_av_qual => [:alt,:ref],
     :t_strand_counts => [:alt_forward,:alt_reverse,:ref_forward,:ref_reverse],
  }
  def initialize file
    @muts = Hash[File.foreach(file).map do |l|
       mut = SomaticIndelOut::Indel.new l.chomp.split(/\t/)
       [ mut.key, mut ]
    end]
  end

  class Indel
    def initialize arr
      @mut = {
        :chrom => arr.shift,
        :start => arr.shift,
        :stop => arr.shift,
        :indel => arr.shift
      }
      read_fields arr
    end

    def read_fields arr
      arr.each do |blob|
        r = blob.match(FIELDS)
        name = r[:field_name].downcase.to_sym
        @mut[ r[:field_name].downcase.to_sym ] = case
        when r[:fields] && (FIELD_KEYS[name] || r[:params])
          Hash[
            (FIELD_KEYS[name] || r[:params].split(%r!/!)).zip(
              r[:fields].split(%r!/!)
            )]
        when r[:fields] && !FIELD_KEYS[name] && !r[:params]
          r[:fields]
        else
          true
        end
      end
    end

    def key
      [ @mut[:chrom], @mut[:start].to_i, @mut[:stop].to_i ]
    end

    def method_missing meth, *args, &block
      @mut[meth] || super(meth, *args, &block)
    end

    def normal_depth
      @normal_depth ||= n_obs_counts[:total].to_i
    end

    def normal_alt_count
      @normal_alt ||= n_obs_counts[:alt].to_i
    end

    def normal_nonref_count
      @normal_nonref_count ||= n_obs_counts[:any].to_i
    end

    def normal_allelic_depth
      [ normal_depth - normal_nonref_count, normal_alt_count ].join ","
    end

    def tumor_depth
      @tumor_depth ||= t_obs_counts[:total].to_i
    end

    def tumor_alt_count
      @tumor_alt_count ||= t_obs_counts[:alt].to_i
    end

    def tumor_nonref_count
      @tumor_nonref_count ||= t_obs_counts[:any].to_i
    end

    def tumor_allelic_depth
      [ tumor_depth - tumor_nonref_count, tumor_alt_count ].join ","
    end
  end

  def [] ind
    @muts[ind]
  end

  def each
    @muts.each do |loc,mut|
      yield mut
    end
  end
end
