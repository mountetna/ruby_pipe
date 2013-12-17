require 'pipeline'
require 'hash_table'

module Utility
  class Config
    extend Pipeline::Config
    include Pipeline::BaseConfig

    def_var :feature_name do cohort_name end

    dir_tree({
      ":output_dir" => {
        "@sample_name.barcode_count.txt" => :barcode_count
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
          genes.each do |gene,barcodes|
            f.puts "#{gene}\t#{barcodes.map{|b,c| "#{b}:#{c}"}.join("\t")}"
          end
        end
      end
    end
  end
  class Cellecta
    include Pipeline::Script
    runs_steps :decipher_barcodes
    def_module :default, :decipher_barcodes => true
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
