require 'pipeline'
require 'hash_table'
require 'fastq'

class CellectaModule < HashTable
  class CellectaBarcode < HashTable::Row
    attr_reader :count
    def add_read
      @count ||= 0
      @count += 1
    end
  end
  index :hugo_symbol

  def initialize opts = {}
    super opts
    @unknown = 0
  end

  def barcode_alias
    @barcode_alias ||= build_barcode_alias
  end

  def build_barcode_alias
    bca = {}
    each do |m|
      m.barcode.size.times do |i|
        x = m.barcode.dup
        x[i] = 'N'
        bca[x] = m
      end
      bca[m.barcode] = m
    end
    bca
  end

  def parse_barcodes(fastq_file)
    fastq = Fastq.new fastq_file
    fastq.each_read do |read|
      barcode = read.seq[0..17].reverse.tr("ATGC", "TACG")
      m = find_barcode(barcode)
      if !m
        unknown_barcode barcode
      else
        m.add_read
      end
    end
  end

  def unknown_barcode barcode
    @unknown += 1
  end

  def find_barcode barcode
    barcode_alias[barcode]

    # find_by_regexp barcode
  end

  def find_by_regexp barcode
    barcode_pattern = Regexp.new barcode.gsub(/N/,'.')
    # scan for it
    matching = select do |m|
      m.barcode =~ barcode_pattern
    end
    matching.first if matching.size == 1
  end

  def write_barcode_count file
    File.open(file,"w") do |f|
      f.puts "gene\tbarcode_id\tcount"
      f.puts "unknown\t?\t#{@unknown}"
      index[:hugo_symbol].entries.each do |gene|
        index[:hugo_symbol][gene].each do |bc|
          f.puts "#{gene}\t#{bc.barcodeid}\t#{bc.count}"
        end
      end
    end
  end
end

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
        #s.extend_with :genesets => make_genesets(s)
        s.replicates.each do |r|
          r.add_member :replicate_name, "r#{r.index}" if !r.replicate_name
        end
      end
    end

    # def_var :feature_name do cohort_name end
    def_var :fdr_cutoff do 1; end
    def_var :test_samples do samples.select{|s| s.test_against} end
    def_var :permutation_fragments do genesets.map{|gs| permutation_fragment(gs)} end

    def_var :diff_exps do samples.collect(&:diff_exp).compact.flatten end
    def_var :sample_replicate_name do |r| "#{(r || job_item).property :sample_name}.#{(r || job_item).property :replicate_name}" end

    dir_tree({
      ":output_dir" => {
        "@sample_name.permutation_result.txt" => :permutation_result,
        "@sample_name.@normal_name.diff_exp" => :diff_exp_table,
        "pdfs" => {
          "@sample_name.high_pval_plot.pdf" => :high_pval_plot,
          "@sample_name.low_pval_plot.pdf" => :low_pval_plot,
          "@sample_name.log_fold_change_plot.pdf" => :log_fold_change_plot
        }
      },
      ":scratch_dir" => {
        "@cohort_name" => {
          "@cohort_name.barcode_count.txt" => :coverage_table,
        },
        "@sample_name" => {
          "@sample_name.gene_list.txt" => :gene_list,
          "@sample_name.@replicate_name.barcode_count.txt" => :barcode_count,
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
    runs_on :samples, :replicates
    audit_report :cellecta_module, :inputs

    class DecipherBarcode
      include Pipeline::Task
      requires_file :inputs, :cellecta_module
      outs_file :barcode_count


      def run
        mod = CellectaModule.new.parse(config.cellecta_module)

        config.inputs.each do |gz|
          mod.parse_barcodes gz
        end
        mod.write_barcode_count config.barcode_count
      end
    end
  end

  class MakeCoverageTable
    include Pipeline::Step
    runs_tasks :coverage_table
    runs_on :cohort

    class CoverageTable
      include Pipeline::Task
      requires_files :samples__replicates__barcode_counts
      outs_file :coverage_table
      def run
        names = config.samples.map do |s|
          s.replicates.map do |r|
            config.sample_replicate_name(r).to_sym
          end
        end.flatten
        combined = HashTable.new :columns => [ :gene, :barcode_id, names ].flatten, :index => :barcode_id
        config.samples.each do |s|
          s.replicates.each do |rep|
            name = config.sample_replicate_name(rep).to_sym
            counts_file = HashTable.new.parse config.barcode_count(rep)
            counts_file.each do |bc|
              line = combined.index[:barcode_id][bc.barcode_id]
              if line.count > 0
                line.first[ name ] = bc.count if bc.count
              else
                new_bc = { :barcode_id => bc.barcode_id, :gene => bc.gene }
                new_bc.update Hash[names.zip [0]*names.size]
                new_bc[name] = bc.count if bc.count
                combined << new_bc
              end
            end
          end
        end
        combined.print config.coverage_table
      end
    end
  end

  class DeseqBarcodes
    include Pipeline::Step
    runs_on :diff_exps
    runs_tasks :deseq

    class Deseq
      include Pipeline::Task
      requires_file :coverage_table
      outs_file :diff_exp_table
      def run
        r_script :deseq, :doDeseq, config.coverage_table, config.sample_name, config.normal_name, config.fdr_cutoff, config.diff_exp_table or error_exit "Could not run DESeq"
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
    runs_steps :decipher_barcodes, :make_coverage_table, :deseq_barcodes, :preprocess_barcodes, :permutation_test, :assemble_fragments, :graph_cellecta_results
    def_module :default, {
      :decipher_barcodes => true,
      :make_coverage_table => true,
      :deseq_barcodes => true,
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
