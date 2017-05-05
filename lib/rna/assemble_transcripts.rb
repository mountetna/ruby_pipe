#!/usr/bin/env ruby
require 'tempfile'
require 'set'

module Rna
  class AssembleTranscripts
    include Pipeline::Step
    runs_tasks :make_fpkm_table, :make_coverage_table

    class MakeFpkmTable
      include Pipeline::Task
      requires_files :gene_trackings
      outs_file :fpkm_table

      def run
        summary = HashTable.new columns: [ :gene_id ] + config.samples__replicatess.map{ |rep| config.sample_replicate_name(rep).to_sym }, index: [ :gene_id ]
        config.samples__replicatess.each do |rep|
          gene_exps = HashTable.new.parse config.gene_tracking(rep)
          name = config.sample_replicate_name(rep).to_sym
          gene_exps.each do |exp|
            if summary.index[:gene_id][exp.gid].count == 0
              summary << { 
                gene_id: exp.gene_id, 
                name => gene.fpkm
              }
            else
              entry = summary.index[:gene_id][exp.gene_id].first
              entry.update name => exp.fpkm
            end
          end
        end
        summary.print config.fpkm_table
      end
    end

    class MakeCoverageTable
      include Pipeline::Task
      requires_files :transcripts_covs
      outs_file :coverage_table

      def run
        combined = {}
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            genes = HashTable.new config.transcripts_cov(rep), :header => [ :gene_id, :coverage ]
            skip = [ "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique" ]
            genes.each do |l|
              next if skip.include? l.gene_id
              combined[l.gene_id] ||= {}
              combined[l.gene_id][config.sample_replicate_name(rep)] = l.coverage
            end
          end
        end
        File.open config.coverage_table, "w" do |f|
          f.puts "gene_id\t#{config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r) } }.flatten.join("\t")}"
          combined.each do |gid,g|
            f.puts "#{gid}\t#{config.samples.map{|s| s.replicates.map{|r| g[config.sample_replicate_name(r)] } }.flatten.join("\t")}"
          end
        end
      end
    end
  end
  class AssembleRsemTranscripts
    include Pipeline::Step
    runs_tasks :make_tpm_table, :make_coverage_table, :make_gene_exp_table, :make_rna_seq_table
    resources memory: "20gb"

    class MakeTpmTable
      include Pipeline::Task
      requires_files :samples__replicates__rsem_genes_resultss
      outs_file :tpm_table

      def run
        names = config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r).to_sym} }.flatten
        summary = HashTable.new columns: [ :gene_id ] + names, index: [ :gene_id ]
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            name = config.sample_replicate_name(rep).to_sym
            log_info "Parsing replicate #{name}"
            gene_exps = HashTable.new.parse config.rsem_genes_results(rep)
            gene_exps.each do |exp|
              if summary.index[:gene_id][exp.gene_id].count == 0
                summary << { 
                  gene_id: exp.gene_id, 
                  name => exp.tpm
                }
              else
                entry = summary.index[:gene_id][exp.gene_id].first
                entry.update name => exp.tpm
              end
            end
          end
        end
        summary.print config.tpm_table
      end
    end

    class MakeCoverageTable
      include Pipeline::Task
      requires_files :samples__replicates__rsem_genes_resultss
      outs_file :coverage_table

      def run
        names = config.samples.map{|s| s.replicates.map{|r| config.sample_replicate_name(r).to_sym} }.flatten
        summary = HashTable.new columns: [ :gene_id ] + names, index: [ :gene_id ]
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            gene_exps = HashTable.new.parse config.rsem_genes_results(rep)
            name = config.sample_replicate_name(rep).to_sym
            gene_exps.each do |exp|
              if summary.index[:gene_id][exp.gene_id].count == 0
                summary << { 
                  gene_id: exp.gene_id, 
                  name => exp.expected_count
                }
              else
                entry = summary.index[:gene_id][exp.gene_id].first
                entry.update name => exp.expected_count
              end
            end
          end
        end
        summary.print config.coverage_table
      end
    end

    class MakeGeneExpTable
      include Pipeline::Task
      requires_files :samples__replicates__rsem_genes_resultss
      outs_file :gene_exp_table

      def hugo_name ensembl_id
        hgnc_entry = hugo_table.index[:ensembl_gene_id][ ensembl_id ].first
        hgnc_entry && hgnc_entry.status != "Entry Withdrawn" ? hgnc_entry.symbol : nil
      end

      def hugo_table
        @hugo_table ||= HashTable.new(index: [ :ensembl_gene_id ]).parse(config.hgnc_complete)
      end

      def run
        all_gene_exp = HashTable.new columns: [ :rna_seq, :ensembl_name, :hugo_name, :read_counts, :expression ]
        config.samples.each do |sample|
          sample.replicates.each do |rep|
            gene_exps = HashTable.new.parse config.rsem_genes_results(rep)
            gene_exps.each do |exp|
              all_gene_exp << { 
                rna_seq: config.sample_name(sample),
                ensembl_name: exp.gene_id,
                hugo_name: hugo_name(exp.gene_id),
                read_counts: exp.expected_count,
                expression: exp.tpm
              }
            end
          end
        end
        all_gene_exp.print config.gene_exp_table
      end
    end

    class MakeRnaSeqTable
      include Pipeline::Task
      requires_files :samples__replicates__rsem_genes_resultss
      outs_file :rna_seq_table

      LOG2_EISENBERG_CUTOFF = {  
        "ENSG00000143612" =>    2.369688,
        "ENSG00000130724" =>    3.496396,
        "ENSG00000105220" =>    2.061996,
        "ENSG00000126067" =>    2.124033,
        "ENSG00000159377" =>    4.308065,
        "ENSG00000075785" =>    8.398825,
        "ENSG00000129625" =>    1.496904,
        "ENSG00000100028" =>    0.000000,
        "ENSG00000165280" =>    2.543607,
        "ENSG00000111237" =>    5.288628
      }

      def eisenberg_score gene_exp
        LOG2_EISENBERG_CUTOFF.count do |gene_id, cutoff|
          Math.log(gene_exp.index[:gene_id][gene_id].first.tpm+0.1,2) > cutoff
        end
      end

      def run
        all_rna_seq = HashTable.new columns: [ :tube_name, :sample,
                                               :rna_seq_plate,
                                               :expression_type,
                                               :transcriptome_build,
                                               :read_count, :duplicate_count,
                                               :mapped_count,
                                               :intergenic_count,
                                               :introns_count, :utr_count,
                                               :coding_count, :mt_coding_count,
                                               :rrna_count, :mt_rrna_count,
                                               :median_3prime_bias,
                                               :median_cv_coverage,
                                               :eisenberg_score ]

        ensembl_pc_genes = HashTable.new.parse(config.ensembl_pc_genes)
        non_mt_genes = Set.new(ensembl_pc_genes.reject { |gene| gene.chrom == "MT" }.map(&:gene_id))
        mt_genes = Set.new(ensembl_pc_genes.select { |gene| gene.chrom == "MT" }.map(&:gene_id))

        config.samples.each do |sample|
          sample.replicates.each do |rep|
            puts config.sample_name(sample)
            tube_name = config.sample_name(sample)
            actual_sample_name = tube_name =~ /.rna./ ? tube_name.sub(/.rna.*/,'') : nil
            flags = Flagstat.new(config.qc_flag(rep))
            rnaseq_metrics = PicardMetrics.new.parse(config.qc_rnaseq(rep))
            bwa_rnaseq_metrics = PicardMetrics.new.parse(config.qc_rnaseq(rep))
            gene_exps = HashTable.new(index: [ :gene_id ], types: { expected_count: :float, tpm: :float }).parse config.rsem_genes_results(rep)
            rrna_metrics = HashTable.new( index: [:stat], columns: [ :stat, :count ], types: { count: :int }, parse_mode: :noheader).parse(config.qc_rrna_metrics(rep))
            all_rna_seq << {
              tube_name: tube_name,
              sample: actual_sample_name,
              rna_seq_plate: config.cohort_name.capitalize,
              expression_type: "TPM",
              transcriptome_build: "Ensembl GRCh38.85",
              read_count: flags.total,
              duplicate_count: flags.duplicates,
              mapped_count: flags.mapped,
              intergenic_count: (rnaseq_metrics.sections[:rna_seq_metrics].intergenic_bases + bwa_rnaseq_metrics.sections[:rna_seq_metrics].intergenic_bases)/config.read_size,
              introns_count: (rnaseq_metrics.sections[:rna_seq_metrics].intronic_bases + bwa_rnaseq_metrics.sections[:rna_seq_metrics].intronic_bases)/config.read_size,
              utr_count: rnaseq_metrics.sections[:rna_seq_metrics].utr_bases/config.read_size,

              coding_count: gene_exps.select{|gexp| non_mt_genes.include?(gexp.gene_id)}.map(&:expected_count).inject(&:+),
              mt_coding_count: gene_exps.select{|gexp| mt_genes.include?(gexp.gene_id)}.map(&:expected_count).inject(&:+),

              rrna_count: rrna_metrics.index[:stat]["rRNA_count"].first.count,
              mt_rrna_count: rrna_metrics.index[:stat]["mt_rRNA_count"].first.count,
              median_3prime_bias: rnaseq_metrics.sections[:rna_seq_metrics].median_3prime_bias,
              median_cv_coverage: rnaseq_metrics.sections[:rna_seq_metrics].median_cv_coverage,

              eisenberg_score: eisenberg_score(gene_exps)
            }
          end
        end
        all_rna_seq.print config.rna_seq_table
      end
    end
  end
end
