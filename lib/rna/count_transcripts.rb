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
  class RsemCount
    include Pipeline::Step
    runs_tasks :rsem_count
    has_tasks :rsem_single_count
    audit_report :sample_replicate_name
    resources :threads => 12
    runs_on :samples, :replicates

    class RsemCount
      include Pipeline::Task
      requires_file :input_fastq1s, :input_fastq2s
      outs_file :output_bam, :output_bai, :rsem_genes_results, :rsem_isoforms_results
      
      def run
        args = {
          fq1s: rsem_format(*config.input_fastq1s),
          fq2s: rsem_format(*config.input_fastq2s),
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
        FileUtils.mv config.rsem_scratch_genome_sorted_bam, config.output_bam or error_exit "Could not move file"
        FileUtils.mv config.rsem_scratch_genome_bai, config.output_bai or error_exit "Could not move file"

        # Move gene and isoform counts to output
        FileUtils.mv config.rsem_scratch_genes_results, config.rsem_genes_results or error_exit "Could not move file"
        FileUtils.mv config.rsem_scratch_isoforms_results, config.rsem_isoforms_results or error_exit "Could not move file"

        FileUtils.rm_rf config.rsem_scratch_dir or error_exit "Could not remove RSEM scratch dir"
      end
    end
    class RsemSingleCount
      include Pipeline::Task
      requires_file :input_fastq1s
      outs_file :rsem_scratch_genes_results, :rsem_scratch_isoforms_results, :rsem_scratch_genome_unsorted_bam
      
      def run
        args = {
          fq1s: rsem_format(*config.input_fastq1s),
          reference: rsem_format(config.reference_rsem),
          sample_name: config.sample_replicate_name
        }
        rsem :calculate_expression, 
            temporary_folder: File.realdirpath(config.rsem_tmp_dir),
            num_threads: config.threads,
            output_genome_bam: true,
            args: args, output: config.rsem_scratch_dir or error_exit "Could not run RSEM"
      end
    end
  end
  class BwaAlign
    include Pipeline::Step
    runs_on :samples, :replicates
    runs_tasks :collect_unmapped, :align_unmapped
    resources :threads => 12, memory: "3gb"

    class CollectUnmapped
      include Pipeline::Task
      requires_file :output_bam
      outs_file :unaligned1_fastq_gz, :unaligned2_fastq_gz

      def run
        log_info "Culling unaligned reads"

        picard :view_sam,
          :I => config.output_bam,
          :ALIGNMENT_STATUS => :Unaligned,
          :out => config.unaligned_sam or error_exit "picard view_sam failed"

        picard :sam_to_fastq,
          :I => config.unaligned_sam,
          :FASTQ => config.unaligned1_fastq,
          :SECOND_END_FASTQ => config.unaligned2_fastq or error_exit "picard sam_to_fastq failed"

        run_cmd "gzip -c #{config.unaligned1_fastq} > #{config.unaligned1_fastq_gz}" or error_exit "Could not gzip fastq file"
        run_cmd "gzip -c #{config.unaligned2_fastq} > #{config.unaligned2_fastq_gz}" or error_exit "Could not gzip fastq file"

        FileUtils.rm config.unaligned_sam or error_exit "Could not delete sam file."
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
