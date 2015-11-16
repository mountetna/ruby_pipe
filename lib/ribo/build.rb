module Ribo
  class Build
    include Pipeline::Step
    runs_tasks :create_transcript_model
    runs_on :cohort


    class CreateTranscriptModel
      include Pipeline::Task
      requires_file :reference_gtf
      outs_file :transcript_model_gtf

      # This transcript model has the following features:
      #   This assumes you are using an ensembl-derived GTF
      #   Only expressed, coding genes should be included. 
      # Then, we define the following regions, per transcript:
      #   5' UTR - any reads falling in the 5' UTR up to 20 bp before the
      #            first ATG
      #   Start site - any reads falling in the region 20 bp before the start site
      #   to 45 bp after the start site.
      #   ORF - any reads falling in the region
      #   end site - any reads falling in the region 20bp before the stop site to 45 bp after the stop site
      #   3' UTR - 45bp after the stop site to the end.
      #
      
      
      class ExonException < StandardError
      end
      def get_exon_region( txp, name, start, stop )
        if start >= stop || stop >= txp.exon_pos.length || start >= txp.exon_pos.length
          log_info "Region has improper bounds! #{name} #{txp.name} #{start} #{stop} #{txp.exon_pos.length}" 
        end

        low, high = [ start, stop ].map do |p| txp.exon_pos[p].pos; end.sort
        region = GenomicLocus::Region.new txp.seqname, low, high
        txp.exons.map do |exon|
          next unless exon.overlaps? region
          if region.contains? exon
            exon.clone do |c|
              c.feature = name
            end
          else
            exon.clone do |c|
              c.feature = name
              c.stop = [ region.stop, c.stop ].min
              c.start = [ region.start, c.start ].max
            end
          end
        end.compact
      end

      def run
        log_info "Loading GTF"
        g = GTF.new(index: [ :gene_name ]).parse( config.reference_gtf )

        genes = g.index[:gene_name].entries
        log_info "GTF has #{genes.count} genes."

        # okay, you should be able to do this just with exons + cds fragments.

        unified = g.wrap([])

        genes.each do |name|
          regions = { }
          config.transcript_model_regions.each do |name|
            regions[name] = []
          end
          g.gene(name).transcripts.each do |txp|
            next if txp.cds.count == 0 || !txp.is_ccds?
            puts "#{txp.gene_name} #{txp.transcript_id}"

            start_pos = txp.translation_start_pos
            start_index = txp.exon_pos.index do |l| l.pos == start_pos.pos end

            stop_pos = txp.translation_stop_pos
            stop_index = txp.exon_pos.index do |l| l.pos == stop_pos.pos end
              
            next if stop_index - start_index < 68 # there is no valid orf here for our purposes.

            # okay, your regions are:

            regions[:utr5].concat get_exon_region(txp, :utr5, 0, start_index - 20 )
            regions[:start].concat get_exon_region(txp, :start, start_index - 19, start_index + 45)
            regions[:orf].concat get_exon_region(txp, :orf, start_index + 46, stop_index - 20)

            stop_end = [ stop_index + 45, txp.exon_pos.count - 1].min
            regions[:stop].concat get_exon_region(txp, :stop, stop_index - 19, stop_end )
            regions[:utr3].concat get_exon_region(txp, :utr3, stop_end, txp.exon_pos.count - 1 )
          end
          regions.each do |name,list|
            r = g.wrap(list).sort_by &:start
            r = r.flatten
            unified.concat r
          end
        end

        unified.print config.transcript_model_gtf
      end
    end
  end
end
