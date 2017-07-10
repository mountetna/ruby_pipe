#!/usr/bin/env ruby
require 'hash_table'
require 'sam'
module Rna
  class CufflinksCount
    include Pipeline::Step
    runs_tasks :cufflink, :format_transcript, :sort_sam, :count_coverage
    has_tasks :cufflink, :cufflink_denovo, :format_transcript, :sort_sam, :count_coverage
    audit_report :sample_replicate_name
    resources :threads => 12
    runs_on :samples, :replicates

    class Cufflink
      include Pipeline::Task
      requires_files :replicate_bam
      dumps_file :transcripts_gtf

      def run
        log_info "Running cufflinks"
        cufflinks :GTF => config.reference_gtf, :bam => config.replicate_bam, :output_dir => config.cufflinks_scratch or error_exit "Cufflinks failed."
      end
    end

    class CufflinkDenovo
      include Pipeline::Task
      requires_files :replicate_bam
      dumps_file :transcripts_gtf

      def run
        log_info "Running cufflinks"
        cufflinks :max_bundle_frags => 10_000_000, :bam => config.replicate_bam, :output_dir => config.cufflinks_scratch or error_exit "Cufflinks failed."
      end
    end

    class FormatTranscript
      include Pipeline::Task
      requires_file :transcripts_gtf
      outs_file :output_gtf

      def run
        log_info "Formatting transcript gtf"
        FileUtils.cp config.transcripts_gtf, config.output_gtf
      end
    end

    class SortSam
      include Pipeline::Task
      requires_file :replicate_bam
      dumps_file :coverage_sam

      def run
        run_cmd "samtools view -b -h -q 1 -f 2 #{config.replicate_bam} | samtools sort -on - #{config.coverage_sam.sub(/.sam/,"")} | samtools view -h - > #{config.coverage_sam}" or error_exit "Could not sort bam correctly."
      end
    end

    class CountCoverage
      include Pipeline::Task
      requires_file :coverage_sam
      outs_file :transcripts_cov

      def run
        log_info "Mapping coverage to reference genes"
        htseq_count :input => config.coverage_sam, :gtf => config.reference_gtf, :type => "exon", :out => config.transcripts_cov or error_exit "Computing normal coverage failed."
      end
    end
  end
  class DepleteRibo
    include Pipeline::Step
    runs_tasks :soak_ribo, :cull_non_ribo, :collect_rrna_metrics
    has_tasks :soak_ribo_single_end, :cull_non_ribo_single_end, :soak_ribo, :cull_non_ribo, :collect_rrna_metrics
    runs_on :samples, :replicates
    resources :threads => 12
    resources memory: "50gb"

    class SoakRibo
      include Pipeline::Task
      requires_files :input_fastq1s, :input_fastq2s
      dumps_file :rrna_sam

      def run
        log_info "Pairing aligned reads"
        bwa_mem fq1:  config.input_fastq1s,
          fq2:  config.input_fastq2s,
          index: config.rrna_bwa_idx, 
          min_score: 23, out: config.rrna_sam or error_exit "BWA mem failed"
      end
    end

    class SoakRiboSingleEnd
      include Pipeline::Task
      requires_files :input_fastq1s
      dumps_file :rrna_sam

      def run
        log_info "Pairing aligned reads"
        bwa_mem fq1:  config.input_fastq1s, 
          index: config.rrna_bwa_idx, 
          min_score: 23, 
          out: config.rrna_sam or error_exit "BWA mem failed"
      end
    end

    class CullNonRibo
      include Pipeline::Task
      requires_file :rrna_sam
      dumps_file :non_rrna1_fastq_gz, :non_rrna2_fastq_gz, :rrna_bam

      def run
        log_info "Culling unaligned reads"
        picard :view_sam,
          :I => config.rrna_sam,
          :ALIGNMENT_STATUS => :Unaligned,
          :out => config.non_rrna_sam or error_exit "picard view_sam failed"

        picard :sam_to_fastq,
          :I => config.non_rrna_sam,
          :FASTQ => config.non_rrna1_fastq,
          :SECOND_END_FASTQ => config.non_rrna2_fastq or error_exit "picard sam_to_fastq failed"

        run_cmd "gzip -c #{config.non_rrna1_fastq} > #{config.non_rrna1_fastq_gz}" or error_exit "Could not gzip fastq file"
        run_cmd "gzip -c #{config.non_rrna2_fastq} > #{config.non_rrna2_fastq_gz}" or error_exit "Could not gzip fastq file"
        
        run_cmd  "samtools view -Sb -F 4 #{config.rrna_sam} | samtools sort -o - #{config.rrna_tmp} >#{config.rrna_bam}"
        run_cmd "samtools index #{config.rrna_bam}"

        File.unlink config.non_rrna1_fastq
        File.unlink config.non_rrna2_fastq
        File.unlink config.non_rrna_sam
        File.unlink config.rrna_sam
      end
    end

    class CullNonRiboSingleEnd
      include Pipeline::Task
      requires_file :rrna_sam
      dumps_file :non_rrna1_fastq_gz, :rrna_bam

      def run
        log_info "Culling unaligned reads"
        picard :view_sam,
          :I => config.rrna_sam,
          :ALIGNMENT_STATUS => :Unaligned,
          :out => config.non_rrna_sam or error_exit "picard view_sam failed"

        picard :sam_to_fastq,
          :I => config.non_rrna_sam,
          :FASTQ => config.non_rrna1_fastq or error_exit "picard sam_to_fastq failed"

        run_cmd "gzip -c #{config.non_rrna1_fastq} > #{config.non_rrna1_fastq_gz}" or error_exit "Could not gzip fastq file"
        
        run_cmd  "samtools view -Sb -F 4 #{config.rrna_sam} | samtools sort -o - #{config.rrna_tmp} >#{config.rrna_bam}"
        run_cmd "samtools index #{config.rrna_bam}"

        File.unlink config.non_rrna1_fastq
        File.unlink config.non_rrna_sam
        File.unlink config.rrna_sam
      end
    end

    class CollectRrnaMetrics
      include Pipeline::Task
      requires_file :rrna_bam
      outs_file :qc_rrna_metrics

      def run
        rrna_count = %x{ samtools view -F 4 -c #{config.rrna_bam} }.to_i
        mt_rrna_count = %x{ samtools view -c #{config.rrna_bam} #{config.mt_rna_contigs} }.to_i

        File.open config.qc_rrna_metrics, "w" do |f|
          f.puts "rRNA_count\t#{rrna_count - mt_rrna_count}"
          f.puts "mt_rRNA_count\t#{mt_rrna_count}"
        end

        File.unlink config.rrna_bam
      end
    end
  end
  class RsemCount
    include Pipeline::Step
    runs_tasks :rsem_count, :rsem_mark_duplicates
    has_tasks :rsem_single_count
    audit_report :sample_replicate_name
    resources :threads => 12
    resources memory: "50gb"
    runs_on :samples, :replicates

    class RsemCount
      include Pipeline::Task
      requires_file :non_rrna1_fastq_gz, :non_rrna2_fastq_gz
      dumps_file :rsem_genome_sorted_bam
      outs_file :rsem_genes_results, :rsem_isoforms_results
      
      def run
        ensure_dir config.rsem_scratch_dir
        ensure_dir config.rsem_tmp_dir
        args = {
          fq1s: rsem_format(config.non_rrna1_fastq_gz),
          fq2s: rsem_format(config.non_rrna2_fastq_gz),
          reference: rsem_format(config.reference_rsem),
          sample_name: config.sample_replicate_name
        }
        rsem :calculate_expression, 
          paired_end: true,
          bowtie2: true,
          temporary_folder: File.realdirpath(config.rsem_tmp_dir),
          num_threads: config.threads,
          output_genome_bam: true,
          sort_bam_by_coordinate: true,
          args: args, output: config.rsem_scratch_dir or error_exit "Could not run RSEM"

        # Move bam files to output
        FileUtils.mv config.rsem_scratch_genome_sorted_bam, config.rsem_genome_sorted_bam or error_exit "Could not move genome_sorted_bam"

        # Move gene and isoform counts to output
        FileUtils.mv config.rsem_scratch_genes_results, config.rsem_genes_results or error_exit "Could not move genes_results"
        FileUtils.mv config.rsem_scratch_isoforms_results, config.rsem_isoforms_results or error_exit "Could not move isoforms_results"

        FileUtils.rm_rf config.rsem_scratch_dir or error_exit "Could not remove RSEM scratch dir"

        File.unlink config.non_rrna1_fastq_gz
        File.unlink config.non_rrna2_fastq_gz
      end
    end
    class RsemMarkDuplicates
      include Pipeline::Task
      requires_file :rsem_genome_sorted_bam
      outs_file :output_bam, :output_bai, :duplication_metrics

      def run
        log_info "Mark duplicates"
        opts = {
          :INPUT => config.rsem_genome_sorted_bam,
          :OUTPUT => config.output_bam,
          :METRICS_FILE => config.duplication_metrics, 
          :CREATE_INDEX => :true,
          :REMOVE_DUPLICATES => :false 
        }
        picard :mark_duplicates, opts or error_exit "Mark duplicates failed"
      end
    end
    class RsemSingleCount
      include Pipeline::Task
      requires_file :non_rrna1_fastq_gz
      dumps_file :rsem_genome_sorted_bam
      outs_file :rsem_genes_results, :rsem_isoforms_results

      def run
        ensure_dir config.rsem_scratch_dir
        ensure_dir config.rsem_tmp_dir
        args = {
          fq1s: rsem_format(config.non_rrna1_fastq_gz),
          reference: rsem_format(config.reference_rsem),
          sample_name: config.sample_replicate_name
        }
        rsem :calculate_expression, 
            bowtie2: true,
            temporary_folder: File.realdirpath(config.rsem_tmp_dir),
            num_threads: config.threads,
            output_genome_bam: true,
            sort_bam_by_coordinate: true,
            args: args, output: config.rsem_scratch_dir or error_exit "Could not run RSEM"

        # Move bam files to output
        FileUtils.mv config.rsem_scratch_genome_sorted_bam, config.rsem_genome_sorted_bam or error_exit "Could not move genome_sorted_bam"

        # Move gene and isoform counts to output
        FileUtils.mv config.rsem_scratch_genes_results, config.rsem_genes_results or error_exit "Could not move genes_results"
        FileUtils.mv config.rsem_scratch_isoforms_results, config.rsem_isoforms_results or error_exit "Could not move isoforms_results"

        FileUtils.rm_rf config.rsem_scratch_dir or error_exit "Could not remove RSEM scratch dir"

        File.unlink config.non_rrna1_fastq_gz
      end
    end
  end
  class BwaAlign
    include Pipeline::Step
    runs_on :samples, :replicates
    runs_tasks :collect_unmapped, :align_unmapped
    resources :threads => 12, memory: "50gb"

    class CollectUnmapped
      include Pipeline::Task
      requires_file :output_bam
      outs_file :unaligned1_fastq_gz, :unaligned2_fastq_gz

      def run
        log_info "Culling unaligned reads"

        samtools "view -b -h -f 4", config.output_bam, config.unaligned_bam  or error_exit "Could not extract unmapped reads."

        picard :sam_to_fastq,
          :I => config.unaligned_bam,
          :FASTQ => config.unaligned1_fastq,
          :SECOND_END_FASTQ => config.unaligned2_fastq or error_exit "picard sam_to_fastq failed"

        run_cmd "gzip -c #{config.unaligned1_fastq} > #{config.unaligned1_fastq_gz}" or error_exit "Could not gzip fastq file"
        run_cmd "gzip -c #{config.unaligned2_fastq} > #{config.unaligned2_fastq_gz}" or error_exit "Could not gzip fastq file"

        FileUtils.rm config.unaligned_bam or error_exit "Could not delete sam file."
      end
    end
    class AlignUnmapped
      include Pipeline::Task
      requires_file :unaligned1_fastq_gz, :unaligned2_fastq_gz
      outs_file :bwa_aligned_bam

      def run
        bwa_mem fq1: config.unaligned1_fastq_gz, fq2: config.unaligned2_fastq_gz, out: config.bwa_aligned_sam or error_exit "Could not run BWA"

        picard :sort_sam,
          :I => config.bwa_aligned_sam,
          :OUTPUT => config.bwa_aligned_bam,
          :CREATE_INDEX => :true,
          :SORT_ORDER => :coordinate or error_exit "picard sort_sam failed"

        FileUtils.rm config.bwa_aligned_sam or error_exit "Could not rm bwa sam"

        FileUtils.rm config.unaligned1_fastq_gz or error_exit "Could not rm fastq"
        FileUtils.rm config.unaligned2_fastq_gz or error_exit "Could not rm fastq"
      end
    end
  end
end
