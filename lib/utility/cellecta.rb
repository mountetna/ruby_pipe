require 'pipeline'
require 'hash_table'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def_var :feature_name do cohort_name end
    # def_var :t0 sample do ...
    # def_var :control_sample do ...
    # def_var :experimental_samples do ...

    def_var :test_samples do samples.select{|s| s.test_against} end

    dir_tree({
      ":output_dir" => {
        "@sample_name.barcode_count.txt" => :barcode_count,
        "@sample_name.normalized_counts.txt" => :normalized_counts
      }
      ":scratch_dir" => {
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
      outs_file :barcode_count

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

  # class PermutationTest
  #   include Pipeline::Step
  #   runs_tasks :permuatation_test
  #   runs_on :test_samples_normalized
  #   audit_report :sample_name

  #   class Permutation
  #     include Pipeline::Task
  #     requires_file :normalized_counts
  #     outs_file :fragments_idk

  #     def run
  #     end
  #   end
  # end

  class Cellecta
    include Pipeline::Script
    runs_steps :decipher_barcodes, :preprocess_barcodes
    def_module :default, {
      :decipher_barcodes => true,
      :preprocess_barcodes => true,
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
