
This is a set of pipelines that do various kinds of analysis. Some standard tools are defined and available.

Installing the pipeline
===

ruby_pipe expects certain tools to be in the path when it runs. Amongst them are a ruby installation, an R installation, and various bioinformatic tools such as BWA, Picard, etc.

Invoking a pipeline script
===
Individual pipelines may be invoked through a common script interface based on symbolic links. The file 'run_pipe.rb' will run the pipeline &lt;pipe>&lt;script>.rb if it is symlinked to 'pipe_script'. E.g. if you symlink this file to the name exome_paired_align (somewhere in your path), it will run the script exome/paired_align.rb. Running the pipeline command with no arguments will show you a list of available commands:

    $ exome_paired_align
    Commands:
     generate <cohort_name>                    # Generate a new config file for a cohort of samples.
     audit <config_file.yml>                   # Audit the pipeline to see which steps are complete.
     start <config_file.yml> [<step>]          # Start the pipeline at the beginning or at <step>
     run_step <config_file.yml> <step_name>    # Run just the named step
     list_steps <config_file.yml>              # List steps for this pipeline.
     stop <config_file.yml> [please]           # stop the pipeline, optionally waiting for the current step to finish
     clean <config_file.yml> <scratch|output|list> <step|all> [<task>]# Clean up files from a given run
     timer <config_file.yml>                   # Generate table of time to completion for each step.


Creating a config file
===

The first command, 'generate', will generate a config file  (in YAML format). This
defines the environment for your analysis, most importantly the list
of samples, with sample names, normal names, input fastqs, patient
ids, etc., as well as any flags that ought to be set.

'generate &lt;cohort_name>' will drop you into a mini-shell where you can
setup a configuration and print it to a file. 'generate &lt;config_file>'
will drop you into a mini-shell to continue editing a generated file.
Sometimes I find it easier to use the mini-shell, sometimes I just use
it to make a skeleton and create the rest in my favorite editor (vim).

Here is a short sample config file which will set you up for an exome run:

    ---
    :cohort_name: example
    :pipe: exome
    :script: paired_align
    :frag_size: 300
    :output_dir: "./output"
    :metrics_dir: "./metrics"
    :scratch_dir: "./scratch"
    :log_dir: "./log"
    :samples:
    - :sample_name: normal
      :inputs:
        - :fq1: /data/sample_exp/normal_R1_001.fastq.gz
          :fq2: /data/sample_exp/normal_R2_001.fastq.gz
    - :sample_name: tumor
	  :normal_name: normal
	  :inputs:
	  - :fq1: /data/sample_exp/tumor_R1_001.fastq.gz
        :fq2: /data/sample_exp/tumor_R2_001.fastq.gz
	  - :fq1: /data/sample_exp/tumor_R1_002.fastq.gz
    	:fq2: /data/sample_exp/tumor_R2_002.fastq.gz
	:interval_list: "/resources/intervals/SureSelectHumanAllExonV4_UTRs_hg19___TARGETS.ilist"

The generate shell is helpful in putting together a list of fastq
files for a sample - you can specify this as a pattern:
> input tumor /data/sample_exp/tumor_*.fastq.gz
It will attempt to pair files automatically based on the Illumina
filename formatting; this usually works pretty well. If it doesn't,
setting it up by hand is not terribly tedious.

Modules
===

The pipeline is modular, meaning you don't have to run the whole
thing. For example, if I already have a set of input BAM files, and I
just want to run mutation calling on it, I can do that. I just have to
specify the inputs in my config file:
    :samples:
	- :sample_name: tumor
	  :input_bam: /data/sample_exp/tumor.bam

Then I specify the set of modules I want to run.
:modules: [ find_mutations ]

You can also change the default pipeline behavior this way. For
example, if I want to use BWA MEM and call indels with somatic indel
detector, I can just load those modules:
:modules: [ default, align_bwa_mem, find_mutations_somatic_indel_detector ]

You can get a list of available modules and the steps they invoke this way:
    $ exome_paired_align list_steps example.exome_paired_align.yml available

'list_steps' without 'available' will show you the current list of
steps the pipeline is set up to run.
    $ exome_paired_align list_steps example.exome_paired_align.yml

Pipeline structure
===

The pipeline is broken down on four levels:

-  A **script** is a set of steps - in this case, exome_paired_align is our
script. The actual list of steps is configurable, as we'll see later.

-  A **step** is the basic unit of execution - each step runs as one arrayed
scheduler job. Each scheduler job item for a step is called a trial.
For example, there might be an 'align' step in the pipeline, which
runs one trial for each pair of input fastqs listed in the config
file.

- A **trial** is a single scheduler job item. It runs on a single node and will complete a series of tasks, exiting as soon as one of the tasks fails. If a task has already completed (its output
files already exist) it will be skipped. If a task is missing its
inputs, it will fail and generate an error. Otherwise, the task will
attempt to execute.

- A **task** is a single unitary operation. Ideally it takes a certain input set of files and produces only a certain output set of files. In practice it also produces scratch files which you may or may not want around. An example task might be a BWA alignment - invoke bwa on an input fastq and produce an alignment.

If ANY task in ANY trial generates an error, once the step has
finished (i.e. all of the threads in the job have exited, cleanly or
not), the pipeline will halt and throw an error to the logs.
Otherwise, it will schedule and execute the next step.

Auditing the pipeline
===

Once you have a valid config file, you should be able to run the
'audit' command. This is a useful way to track the current progress of
the pipeline. In this case, we would do:
$ exome_paired_align audit example.exome_paired_align.yml

This will audit the full execution of the pipeline from start to
finish. You might want to pipe to less! (If you do, export LESS=-R to
fix the colors).

The steps are listed in order from first to last, along with the
progress of each trial and each task, with a summary of what tasks are
complete, which files have been made, which are missing, etc. If you
audit initially you should see a lot of red - this means nothing has
been generated yet. Only the first task of the first step should
report that it is ready to run. We'll soon fix that.

The audit is a useful way to debug. You can restrict the audit to a
single step, e.g.:
    $ exome_paired align audit example.exome_paired_align.yml align
if you only want to audit the align step. Usually between this and the
log output (which hopefully includes stack traces of your error) you
should be able to understand what went wrong.

Running the pipeline
===

There are three useful ways to start the pipeline:
    $ exome_paired_align start example.exome_paired_align.yml
This is the basic way, and will start the pipeline at the first step
and run through it.

If things go south, and you need to restart the pipeline from the
middle, no problem:
    $ exome_paired_align start example.exome_paired_align.yml realign
Just give the name of the step you'd like to resume with, and the
pipeline will continue on its merry way. Make sure you audit before
you do this, to make sure that things will actually run.

Finally, if you're really in a rut, and you just want to run this one
step without continuing on to the next one (e.g. if you want to verify
the output of this step before continuing), you can tell the pipeline
to do so:
    $ exome_paired_align run_step example.exome_paired_align.yml realign

There is a stop command:
    $ exome_paired_align stop example.exome_paired_align.yml

This has varying degrees of success, depending on how moody your
scheduler is. Give it a shot; if it doesn't work, you'll have to kill
the jobs by hand (with qdel or canceljob or whatever).

Cleaning up
===

The pipeline generates a LOT of intermediate files, and it can quickly
take up a lot of space if you don't mop up. There is a tool to do some
cleanup, but it's embarrassingly bad, so it might not clean up
everything as well as you'd like. It's usually a good idea to wait
until the end to clean up, mostly because if things mess up, you may
want your temp files so you can resume from the middle.

The first thing you want to do is this:
    $ exome_paired_align clean example.exome_paired_align.yml list <step>

This will produce a summary of how many files have been created by the
step, split into 'scratch' (temp) and 'output' (permanent) categories.
N.B. this division is somewhat arbitrary and needs a lot of cleanup.

You can clean scratch or output files separately:
    $ exome_paired_align clean example.exome_paired_align.yml scratch <step>
The separation is intentional, to prevent accidents. It's also
possible to say 'all' instead of cleaning step-by-step, but I'm
usually not that brave.

Reading logs
===

The logs are copious and may run into the hundreds of megabytes. They
are also in color (I have a fetish), so export LESS=-R to read them.
Usually they write to the ./log directory and have the name format
'&lt;pipe>.&lt;cohort_name>.&lt;step>.&lt;trial>log', e.g.,
"exome.example.mut_det.22.log". ERROR files, which show up in the
working directory, will usually tell you which logs to look at. There
is also a summary log file for the script, e.g.  "exome.example.paired_align.log"
