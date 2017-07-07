# Previously, this pipeline was divided into "Step" and
# "Task", where each "Step" ran a series of tasks across an
# array of job items.
# 
# Instead, let us imagine a pipeline as a series of
# dependencies.  A task is run across a series of input
# objects.
# 
# Typically, a task requires:
#   - A set of input symbols or values A list of input file
#   - locations A set of tools A set of system resource file
#   - locations
# 
# The main question is: when does the task execute?
#   - When it has all of its required inputs.  When its
#   - output is incomplete.
# 
# That is all. There is one job which polls, whose
# responsibility is dispatch. The other jobs keep track of
# their status and log the results to a database.
# 
# The dispatcher has one function: it runs tasks until
# everything refuses to run. Each time it runs, it will
# evaluate which tasks need to run and set up a scheduler
# job to run them.
# 
# The job of the SCRIPT, then is to define the dependency
# relationship of the tasks. This is primarily about
# following intermediate file creation.
# 
# How does this work?
# 
# For example, let us say we want to do mutation detection on
# a bam file.
# 
# The basic process would be something like this:
# 
# For each chromosome in the human genome, scan the BAM for
# mutations.
# For each chromosome, annotate mutations.
# For each sample, assemble chromosomes of somatic mutations.
# 
# in some pseduo-code, here is a simple task.

class Mutect < Pipeline::Task
  input :bam_file, format: :BAM
  input :region
  resource :hg38
  output :mutect_vcf_file
  tool :mutect

  def run
    mutect(etc. bam_file, output: mutect_vcf_file)
  end
end

class MutationDetection < Pipeline::Script
  scratch :chrom_mutect_file, ":scratch_dir/@{sample.name}/@{chrom.name}.mutect.vcf.txt", :VCF
  scratch :mutect_annotated_vcf, ":scratch_dir/@{sample.name}/@{chrom.name}.mutect.annotated.vcf.txt", :VCF
  output :somatic_vcf, ":output_dir/@{sample.name}/@{sample.name}.annotated.vcf"

  across :samples, :chroms do |sample,chrom|
    task :mutect, 
      bam_file: sample.input_bam,
      region: chrom.name,
      mutect_vcf_file: chrom_mutect_file

    task :annotate_vcf,
      raw_vcf: chrom_mutect_file,
      annotated_vcf: mutect_annotated_vcf

    task :filter_vcf,
      raw_vcf: mutect_annotated_vcf,
      filtered_vcf: mutect_filtered_vcf,
      filter: {
      }
  end

  across :samples do |sample|
    task :combine_vcfs,
      vcfs: across(:chroms) { mutect_annotated_vcf },
      vcf: 
  end
end
