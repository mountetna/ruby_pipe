require 'pipeline'
require 'hash_table'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig
    
    def make_genesets s
      if File.exists? gene_list(s)
        genes = File.foreach(gene_list(s)).map &:chomp
        # genes.each_slice(50).map do |gs|
        genes.each_slice(160).with_index.map{ |gs,i| { :genes => gs, :geneset_name => "geneset#{i}" } }
      else
        [ { :genes => [ "DUMMY" ], :geneset_name => "geneset0"} ]
      end
    end

    def init_hook
      #see make_chunks in genome config.rb
      samples.each do |s|
        s.extend_with :genesets => make_genesets(s)
      end
    end

    # def_var :feature_name do cohort_name end
    def_var :test_samples do samples.select{|s| s.test_against} end
    def_var :permutation_fragments do genesets.map{|gs| permutation_fragment(gs)} end

    dir_tree({
      ":output_dir" => {
        "@sample_name.permutation_result.txt" => :permutation_result,
        "pdfs" => {
          "@sample_name.high_pval_plot.pdf" => :high_pval_plot,
          "@sample_name.low_pval_plot.pdf" => :low_pval_plot,
          "@sample_name.log_fold_change_plot.pdf" => :log_fold_change_plot
        }
      },
      ":scratch_dir" => {
        "@sample_name" => {
          "@sample_name.gene_list.txt" => :gene_list,
          "@sample_name.barcode_count.txt" => :barcode_count,
          "@sample_name.normalized_counts.txt" => :normalized_counts,
          "@sample_name.genes_of_interest.txt" => :genes_to_plot,
          "@sample_name.gene_fragments" => {
            "." => :gene_fragment_dir,
            "@sample_name.@geneset_name.genelist_fragment.txt" => :gene_fragment,
            "@sample_name.@geneset_name.permutation_fragment.txt" => :permutation_fragment
          }
        }
      }
    })
  end

  class DecipherBarcodes
    include Pipeline::Step
    runs_tasks :decipher_barcode
    runs_on :samples
    audit_report :cellecta_module, :inputs

    class DecipherBarcode
      include Pipeline::Task
      requires_file :inputs, :cellecta_module
      outs_file :barcode_count, :gene_list

      def get_barcode file
        io = IO.popen("zcat #{file}")
        while read = io.gets
          seq = io.gets
          id = io.gets
          qual = io.gets
          barcode = seq[0..17].reverse.tr("ATGC", "TACG")
          yield barcode
        end
      end

      def build_barcode_alias mod
        bca = {}
        mod.barcode.each do |bc, m|
          bc.size.times do |i|
            x = bc.dup
            x[i] = 'N'
            bca[x] = bc
          end
          bca[bc] = bc
        end
        bca
      end

      def add_barcode genes, bc, mod
        genes[ bc.hugo_symbol ] ||= Hash[mod.hugo_symbol[bc.hugo_symbol].map{|b| [ b.barcodeID, 0 ] }]
        genes[ bc.hugo_symbol ][ bc.barcodeID ] += 1
      end

      def write_gene_list genes
        genes.delete 'unknown'
        File.open(config.gene_list, "w") do |f|
          genes.each do |gene, x|
            f.puts gene
          end
        end
      end

      def run
        mod = HashTable.new config.cellecta_module, :idx => [ :barcode, :hugo_symbol ]
        barcode_alias = build_barcode_alias mod
        genes = { :unknown => { :"?" => 0 } }

        File.open(config.barcode_count,"w") do |f|
          config.inputs.each do |gz|
            get_barcode(gz) do |barcode|
              bc = mod.barcode[barcode_alias[barcode]].first if barcode_alias[barcode]
              if !bc
                genes[:unknown][:"?"] += 1
                next
              end
              add_barcode genes, bc, mod
            end
          end
          f.puts "gene\tbc\tcount"
          genes.each do |gene,barcodes|
            barcodes.each do |bc, count|
              f.puts "#{gene}\t#{bc}\t#{count}"
            end
          end
        end
        write_gene_list genes
      end
    end
  end

  class PreprocessBarcodes
    include Pipeline::Step
    runs_tasks :normalize_counts
    runs_on :test_samples
    audit_report :sample_name

    class NormalizeCounts
      include Pipeline::Task
      requires_file :barcode_count#, :cellecta_module
      outs_file :normalized_counts

      def run
        ctrl_name = config.test_against
        ctrl_samp = config.samples.find{|s| s.sample_name == ctrl_name}
        r_script :cellecta_preprocessing, :normalize_and_log_fold_change, config.normalized_counts, config.barcode_count(ctrl_samp), config.barcode_count
      end
    end
  end

  class PermutationTest
    include Pipeline::Step
    runs_tasks :permutation
    runs_on :test_samples, :genesets
    audit_report :sample_name

    class Permutation
      include Pipeline::Task
      requires_file :normalized_counts
      outs_file :permutation_fragment

      def write_gene_list gene_set, geneset_name
        File.open(config.gene_fragment, "w") do |f|
          gene_set.each do |g|
            if g != 'unknown'
              f.puts g
            end
          end
        end
      end

      def run
        write_gene_list config.genes, config.geneset_name
        r_script :cellecta_permutation_test, :permute, config.normalized_counts, config.gene_fragment, config.permutation_fragment
        File.delete config.gene_fragment
      end
    end
  end

  class AssembleFragments
    include Pipeline::Step
    runs_tasks :glue_fragments
    runs_on :test_samples

    class GlueFragments
      include Pipeline::Task
      requires_file :permutation_fragments
      outs_file :permutation_result

      def run
        open(config.permutation_result, 'a') do |f|
          f.puts "gene\tsecond_highest\thigh_pval\tsecond_lowest\tlow_pval"
          Dir.glob("#{config.gene_fragment_dir}/*") do |fname|
            f.puts File.readlines(fname)[1..-1].join() #removes header line
          end
        end
      end
    end
  end

  class GraphCellectaResults
    include Pipeline::Step
    runs_tasks :graph_pvals, :graph_log_fold_change
    runs_on :test_samples

    class GraphPvals
      include Pipeline::Task
      requires_file :permutation_result
      outs_file :high_pval_plot, :low_pval_plot

      def run
        r_script :cellecta_plots, :make_pval_plots, config.permutation_result, config.sample_name, config.high_pval_plot, config.low_pval_plot
        r_script :cellecta_plots, :make_lfc_plots, config.normalized_counts, config.sample_name, config.log_fold_change_plot
      end
    end

    class GraphLogFoldChange
      include Pipeline::Task
      requires_file :normalized_counts
      outs_file :log_fold_change_plot

      def write_genes_of_interest
        File.open(config.genes_to_plot, "w") do |f|
          config.genes_of_interest.each do |g|
            print "g = #{g}"
            f.puts g
          end
        end
      end

      def run
        write_genes_of_interest
        r_script :cellecta_plots, :make_lfc_plots, config.normalized_counts, config.sample_name, config.log_fold_change_plot, config.genes_to_plot
      end
    end
  end

  class Cellecta
    include Pipeline::Script
    runs_steps :decipher_barcodes, :preprocess_barcodes, :permutation_test, :assemble_fragments, :graph_cellecta_results
    def_module :default, {
      :decipher_barcodes => true,
      :preprocess_barcodes => true,
      :permutation_test => true,
      :assemble_fragments => true,
      :graph_cellecta_results => true
    }

    class ConfigGenerator
      include Pipeline::ConfigGenerator
      include Pipeline::Usage

      def input(args)
        sample = get_sample args.first
        sample[:inputs] = []
        getlines "your fastq files" do |l|
          sample[:inputs] += l.strip.split
        end
      end
      usage "input <sample name>", "Input fastq file for the given sample"

      def module(args)
        config[:cellecta_module] = args.first
      end
      usage "module <file>", "Specify the cellecta module you wish to decipher."
    end
  end
end
