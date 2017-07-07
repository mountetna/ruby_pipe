
Previously, this pipeline was divided into "Step" and
"Task", where each "Step" ran a series of tasks across an
array of job items.

Instead, let us imagine a pipeline as a series of
dependencies.  A task is run across a series of input
objects.

Typically, a task requires:
  - A set of input symbols or values A list of input file
  - locations A set of tools A set of system resource file
  - locations

The main question is: when does the task execute?
  - When it has all of its required inputs.  When its
  - output is incomplete.

That is all. There is one job which polls, whose
responsibility is dispatch. The other jobs keep track of
their status and log the results to a database.

The dispatcher has one function: it runs tasks until
everything refuses to run. Each time it runs, it will
evaluate which tasks need to run and set up a scheduler
job to run them.

The job of the SCRIPT, then is to define the dependency
relationship of the tasks. This is primarily about
following intermediate file creation.

How does this work?

For example, let us say we want to do mutation detection on
a bam file.

The basic process would be something like this:

For each chromosome in the human genome, scan the BAM for
mutations.

For each chromosome, annotate mutations.
For each sample, assemble chromosomes of somatic mutations.

in some pseduo-code, here is a simple task.

class Mutect < Pipeline::Task
  input :bam_file, format: :BAM
  input :region
  resource :hg38
  output :mutect_file
  tool :mutect

  def run
    mutect(etc. bam_file, output: mutect_file)
  end
end

class MutationDetection < Pipeline::Script
  class Sample < Pipeline::Object
    def output_bam
    end
  end
  def script
    samples.each do |sample|
      Mutect.new(bam_file: sample.output_bam,
    end
  end
end
