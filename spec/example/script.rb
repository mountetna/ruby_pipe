require_relative '../pipeline/task'
require_relative '../pipeline/script'

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
  input :region, format: String
  output :mutect_vcf_file, format: :VCF

  tool :mutect
end

class AnnotateVcf < Pipeline::Task
  input :raw_vcf, format: :VCF
  output :annotated_vcf, format: :VCF

  tool :snpeff
end

class FilterVcf < Pipeline::Task
  input :raw_vcf, format: :VCF
  output :filtered_vcf, format: :VCF

  tool :vcftools
end

class CombineVcf < Pipeline::Task
  input :vcfs, format: [ :VCF ]
  output :vcf, format: :VCF

  tool :vcftools
end


# Now a question is - where do samples and chroms come from? They should be objects, of course - but of what sort?
class Sample < Pipeline::Object
  property :input_bam, String
  property :inputs, [ Input ]
end

class MutationDetection < Pipeline::Script
  object :samples, config(:samples)
  object :chroms, config(:genome)

  scratch :chrom_mutect_file, ":scratch_dir/@sample/@chrom.mutect.vcf.txt"
  scratch :mutect_annotated_vcf, ":scratch_dir/@sample/@chrom.mutect.annotated.vcf.txt"

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

  output :somatic_vcf, ":output_dir/@sample/@sample.annotated.vcf"

  across :samples do |sample|
    task :combine_vcfs,
      vcfs: collect(:chroms) { mutect_annotated_vcf },
      vcf: somatic_vcf
  end
end
