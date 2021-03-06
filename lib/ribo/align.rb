module Ribo
  class Align
    include Pipeline::Step
    runs_tasks :clip_fastq, :align_single, :map_reads, :enforce_label
    has_tasks :clip_fastq, :align_single, :map_reads, :enforce_label, :soak_ribo, :cull_non_ribo, :collect_rrna_metrics, :make_nonribo_fastq
    runs_on :fractions
    resources :threads => 6

    class ClipFastq
      include Pipeline::Task
      requires_file :inputs
      dumps_file :clipped_fastq

      def run
        log_info "Clipping adapter from fastq"
        fastx_clipper :adapter => config.adapter,  :out => config.clipped_fastq,
          :in => config.inputs or error_exit "fastx_clipper failed"
      end
    end

    class SoakRibo
      include Pipeline::Task
      requires_files :clipped_fastq
      dumps_file :ribo_sam

      def run
        log_info "Pairing aligned reads"
        bwa_mem fq1:  config.clipped_fastq, index: config.ribo_bwa_idx, min_score: 23, out: config.ribo_sam or error_exit "BWA mem failed"
      end
    end

    class CullNonRibo
      include Pipeline::Task
      requires_file :ribo_sam
      outs_file :non_ribo_sam

      def run
        log_info "Culling unaligned reads"
        picard :view_sam,
          :I => config.ribo_sam,
          :ALIGNMENT_STATUS => :Unaligned,
          :out => config.non_ribo_sam or error_exit "picard view_sam failed"
      end
    end

    class CollectRrnaMetrics
      include Pipeline::Task
      requires_file :ribo_sam, :non_ribo_sam
      dumps_file :qc_rrna_metrics

      def run
        ribo_count = %x{ wc -l #{config.ribo_sam} }.to_i
        non_ribo_count = %x{ wc -l #{config.non_ribo_sam} }.to_i
        header_count = %x{ head -200 #{config.non_ribo_sam} | grep ^@ | wc -l }.to_i

        File.open config.qc_rrna_metrics, "w" do |f|
          f.puts "rRNA_reads\t#{ribo_count - header_count}"
          f.puts "non_rRNA_reads\t#{non_ribo_count - header_count}"
        end
      end
    end
    
    class MakeNonriboFastq
      include Pipeline::Task
      requires_file :non_ribo_sam
      dumps_file :non_ribo_fastq

      def run
        picard :sam_to_fastq,
          :I => config.non_ribo_sam,
          :FASTQ => config.non_ribo_fastq or error_exit "picard sam_to_fastq failed"
      end
    end

    class AlignSingle
      include Pipeline::Task
      requires_file :clipped_fastq
      dumps_file :read_sai

      def run
        log_info "Aligning file with bwa"
        bwa_aln :fq => config.clipped_fastq, :out => config.read_sai or error_exit "bwa aln failed"
      end
    end

    class MapReads
      include Pipeline::Task
      requires_file :read_sai, :clipped_fastq
      dumps_file :mapped_sam

      def run
        log_info "Mapping reads with bwa"
        bwa_single :fq => config.clipped_fastq, :map => config.read_sai, :out => config.mapped_sam or error_exit "bwa samse failed"
      end
    end

    class EnforceLabel 
      include Pipeline::Task
      requires_file :mapped_sam
      outs_file :genome_bam

      def run
        log_info "Enforce read group assignments"
        picard :add_or_replace_read_groups, :INPUT => config.mapped_sam,
                :CREATE_INDEX => :true,
                :SO => :coordinate,
                :OUTPUT => config.genome_bam,
                  :RGID => config.sample_name, :RGLB => config.sample_name,
                  :RGPL => :illumina, :RGPU => :None, :RGSM => config.sample_name,
                  :CN => "Babel@UCSF" or error_exit "Relabel failed"
      end
    end
  end
end
