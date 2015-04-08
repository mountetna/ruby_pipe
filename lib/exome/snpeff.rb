require 'vcf'

class SnpeffVCF < VCF
  attr_reader :tumor_name, :normal_name
  def initialize opts={}
    super opts
    @normal_name = opts[:normal_name]
    @tumor_name = opts[:tumor_name]
  end

  class Line < VCF::Line
    def effects
      @effects ||= parse_effects
    end

    def best_effect
      @best ||= effects.max_by &:score
    end

    def t_alt_count
      @t_alt_count ||= genotype(@table.tumor_name).ad.split(/,/).map(&:to_i)[1]
    end

    def t_ref_count
      @t_ref_count ||= genotype(@table.tumor_name).ad.split(/,/).map(&:to_i)[0]
    end

    def n_alt_count
      @n_alt_count ||= genotype(@table.normal_name).ad.split(/,/).map(&:to_i)[1]
    end

    def n_ref_count
      @n_ref_count ||= genotype(@table.normal_name).ad.split(/,/).map(&:to_i)[0]
    end

    def t_var_freq
      t_alt_count.to_f / (t_alt_count + t_ref_count)
    end

    def n_var_freq
      n_alt_count.to_f / (n_alt_count + n_ref_count)
    end

    private
    class Effect
      ATTRS = [
        :allele,
        :annotation,
        :annotation_impact,
        :gene_name,
        :gene_id,
        :feature_type,
        :feature_id,
        :transcript_biotype,
        :rank,
        :hgvs_c,
        :hgvs_p,
        :cdna,
        :cds,
        :aa,
        :distance,
        :errors 
      ]
      ATTRS.each do |attr|
        attr_reader attr
      end

      attr_accessor :longest

      def initialize values
        values.each do |name,value|
          instance_variable_set("@#{name}", value)
        end
      end

      def cds_length
        cds_array[1] || -1
      end

      def cds_array
        @cds_array ||= cds ? cds.split(/\//).map(&:to_i) : []
      end

      def score
        score=0
        if transcript_biotype=='protein_coding'
          score += 1
        elsif transcript_biotype =~ /_pseudogene/
          score=score-1
        end
        score *= 10

        coding = {
          'missense_variant' => 6, 
          'transcript_amplification' => 6, 
          'splice_donor_variant' => 6, 
          'splice_acceptor_variant' => 6, 
          'initiator_codon_variant' => 5, 
          'frameshift_variant'  =>  5, 
          'stop_gained' => 5, 
          'stop_lost' => 5, 
          'inframe_insertion'  =>  4, 
          'inframe_deletion'  =>  4, 
          'non_conservative_missense_variant' => 4,
          'synonymous_variant' => 4, 
          'transcript_ablation' => 3, 
          'splice_region_variant' => 2, 
          'incomplete_terminal_codon_variant' => 1, 
          'stop_retained_variant' => 1, 
          'coding_sequence_variant' => 1 
        }
        if coding[annotation]
          score += coding[annotation]
          score *= 10
          
          if longest
            score += 1
          end
          score=score*10
        else
          if annotation == 'intron_variant' || annotation=='intergenic_variant'
            score += 1
          elsif annotation =~ /^stream_/
            score -= 1
          end
          score=score*10
        end
        if errors && errors.size > 0
          score -= 1
        end
        if errors =~/WARNING_TRANSCRIPT_NO_START_CODON/
          score -= 100000
        end
        score
      end
    end
    def parse_effects
      return [] unless info[:ANN]

      effects = info[:ANN].split(/,/).map do |txt|
        Effect.new Hash[Effect::ATTRS.zip txt.split(/\|/)]
      end
      max = effects.max_by &:cds_length
      max.longest = true
      effects
    end
  end
end

class MutectSnpeffVCF < SnpeffVCF
  class Line < SnpeffVCF::Line
    def judgement
      info[:JM]
    end

    def covered
      info[:CV]
    end

    def q0_ratio
      info[:MQ0].to_f / total_alt_depth
    end

    def total_alt_depth
      @genotypes.values.inject(0) do |sum,gt|
        sum += gt.ad.split(/,/).last.to_i
        sum
      end
    end
  end
end

class SomaticIndelSnpeffVCF < SnpeffVCF
  class Line < SnpeffVCF::Line
    def alt_count
      @alt_count ||= ad.split(/,/).map(&:to_i)[1]
    end
    def ref_count
      @ref_count ||= ad.split(/,/).map(&:to_i)[0]
    end
    def depth
      @depth ||= dp.to_i
    end
    def alt_mismatch_rate
      nqsmm.to_f
    end
    def alt_mismatch_count
      mm.to_i
    end
    def alt_base_quality; respond_to?(:nqsbq) ? nqsbq.split(/,/)[0].to_f : nil; end
    def alt_map_quality; respond_to?(:mqs) ? mqs.split(/,/)[0].to_f : nil; end
    def alt_freq; alt_count / depth.to_f; end
    def ref_freq; ref_count / depth.to_f; end
    def alt_base_quality; respond_to?(:nqsbq) ? nqsbq.split(/,/)[0].to_f : nil; end
    def alt_map_quality; respond_to?(:mqs) ? mqs.split(/,/)[0].to_f : nil; end
    def alt_mismatch_rate; respond_to?(:nqsmm) ? nqsmm.split(/,/)[0].to_f : nil; end
    def alt_mismatch_count; respond_to?(:mm) ? mm.split(/,/)[0].to_f : nil; end
    def quality; gq.to_i; end
  end
end
