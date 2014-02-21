require 'pipeline'
require 'hash_table'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig
    
    def make_genesets s
      if File.exists? gene_list(s)
        genes = File.foreach(gene_list(s)).map &:chomp
        # genes.each_slice(500).map do |gs|
        genes.each_slice(50).map do |gs|
          { :genes => gs }
        end
      else
        [ { :genes => [ "DUMMY" ]} ]
      end
    end

    def init_hook
      #see make_chunks in genome config.rb
      samples.each do |s|
        s.extend_with :gene_sets => make_genesets(s)
      end
    end

    def_var :feature_name do cohort_name end
    # def_var :t0 sample do ...
    # def_var :control_sample do ...
    # def_var :experimental_samples do ...

    def_var :test_samples do samples.select{|s| s.test_against} end

    dir_tree({
      ":output_dir" => {
        "@sample_name.barcode_count.txt" => :barcode_count,
        "@sample_name.normalized_counts.txt" => :normalized_counts,
        "@sample_name.permutation_results.txt" => :permutation_results
      },
      ":scratch_dir" => {
        "@sample_name.gene_list.txt" => :gene_list,
        "@sample_name" => {
          "@sample_name.gene_fragments" => {
            "." => :gene_fragment_dir,
            "idk" => :dummy_because_this_folder_wont_be_created_without_it
          },
          "@sample_name.output_fragments" => {
            "." => :results_fragment_dir,
            "@sample_name.fragments.txt" => :fragments
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
      requires_file :barcode_count, :cellecta_module
      outs_file :normalized_counts

      def run
        ctrl_name = config.test_against
        ctrl_samp = config.samples.select{|s| s.sample_name == ctrl_name}
        # ctrl_bc_count_loc = ctrl_samp.select{|s| s.barcode_count} #returns nothing
        # ctrl_bc_count_loc = ctrl_samp.collect(&:barcode_count) # throws error
        # print "ctrl name: #{ctrl_name}\n" # works
        # print "ctrl sample: #{ctrl_samp}\n" #works
        # print "ctr; barcode count loc: #{ctrl_bc_count_loc}\n"  
        # print "ctr; barcode count loc: #{ctrl_samp.length}\n" # returns 1
        # print "ctrl samp.barcode_conut: #{ctrl_samp.barcode_count}\n" # throws error 'undefined method barcode_count for array' 
        # print "ctr; barcode count loc[0]: #{ctrl_samp[0]}\n" # returns #<Pipeline::SampleObject:0x0000000227f0b8>
        # print "ctr; barcode count loc[0].barcode_count: #{ctrl_samp[0].barcode_count}\n" #returns blank
        # print "ctr; barcode count loc: #{ctrl_samp[0].select{|s| s.barcode_count}}\n" # returns blank

        ctrl_bc_count_loc = "./output/#{ctrl_name}.barcode_count.txt"

        r_script :cellecta_preprocessing, :normalize_and_log_fold_change, config.normalized_counts, ctrl_bc_count_loc, config.barcode_count
      end
    end
  end

  class PermutationTest
    include Pipeline::Step
    runs_tasks :permutation
    runs_on :test_samples, :gene_sets
    audit_report :sample_name

    class Permutation
      include Pipeline::Task
      requires_file :normalized_counts
      outs_file :fragments, :dummy_because_this_folder_wont_be_created_without_it

      def write_gene_list gene_set, job_index
        fragment_dir = config.gene_fragment_dir[0...-1] #get rid of period at the end...
        File.open("#{fragment_dir}#{job_index}_genes.txt", "w") do |f|
          gene_set.each do |g|
            if g != 'unknown'
              f.puts g
            end
          end
        end
      end

      def run
          # print "config.genes #{config.genes}\n"
          # print "config.results_fragment_dir #{config.results_fragment_dir}\n"
          # print "config.job_index #{config.job_index}\n"
          # print "config.gene_fragment_dir #{config.gene_fragment_dir}\n"
          write_gene_list config.genes, config.job_index
          r_script :cellecta_permutation_test, :permute, config.normalized_counts, config.job_index, config.gene_fragment_dir, config.results_fragment_dir 
      end
    end
  end

  class Cellecta
    include Pipeline::Script
    runs_steps :decipher_barcodes, :preprocess_barcodes, :permutation_test
    def_module :default, {
      :decipher_barcodes => true,
      :preprocess_barcodes => true,
      :permutation_test => true
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
