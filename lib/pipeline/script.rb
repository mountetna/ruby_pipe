require 'fileutils'
module Pipeline
  module Script
    include Pipeline::Base
    include Pipeline::Logger
    include Pipeline::Usage
    include Pipeline::Scheduling

    module ClassMethods 
      attr_reader :steps, :modules
      def runs_steps(*step_list)
        @steps = step_list
      end

      def def_module mod, steps
        @modules ||= {}
        @modules[mod] = steps
      end

      def expand_module m
        modules[m].map do |n,disp|
          modules[n] ? expand_module(n) : { n => disp }
        end.reduce :merge
      end
    end

    def self.included(base)
      base.extend(ClassMethods) if ClassMethods
    end

    def modules
      @modules ||= config.modules.inject({}) do |mods,m|
        mods.update self.class.expand_module m.to_sym
      end
    end

    def steps
      @steps ||= self.class.steps.map{ |s| modules[s] ? s : nil }.compact
    end

    def exclude_task? task
      nil
    end

    def initialize(d=nil)
      @defaults = d
    end

    def config
      @config ||= self.class.sister_class(:config).new(self,@defaults)
    end

    def run_action(*args)
      # run the current action
      send config.action, args
    end

    def schedule(args)
      # handle any previous steps
      curr_step = config.step

      create_step(config.prev_step).complete if config.prev_step

      return if config.single_step

      # schedule the execution for the current step
      job, splits = create_step(curr_step).setup_exec

      # schedule the scheduler for the next step
      create_step(config.next_step).setup_scheduler(job, splits) if config.next_step
    end

    def create_step(s)
      self.class.sister_class(s).new(self, modules[s])
    end

    def exec(args)
      step = create_step(config.step)

      # clean up if you manage to get through the whole step.
      step.cleanup if step.exec 
    end

    def start_pipe(step)
      FileUtils.rm(config.error_file) if File.exists?(config.error_file)
      config.set_config :scheduler, scheduler_type
      job, splits = create_step(step).setup_exec
      create_step(config.next_step).setup_scheduler(job,splits) if config.next_step
    end

    usage "list_steps <config_file.yml>", "List steps for this pipeline."
    usage "generate <cohort_name>", "Generate a new config file for a cohort of samples."
    usage "audit <config_file.yml>", "Audit the pipeline to see which steps are complete."
    usage "start <config_file.yml> [<step>]", "Start the pipeline at the beginning or at <step>"
    usage "run_step <config_file.yml> <step_name>", "Run just the named step"
    usage "stop <config_file.yml> [please]", "stop the pipeline, optionally waiting for the current step to finish"
    usage "clean <config_file.yml> <scratch|output|list> <step|all> [<task>]", "Clean up files from a given run"

    def generate(args)
      self.class.daughter_class(:config_generator).new *args
    end

    def set_config args, act = nil
      ENV['CONFIG'] = args.shift
      abort "Config file not found.".red.bold if !File.exists? ENV['CONFIG']
      config.set_config :action, (act || :init)
      args
    end

    def run_step(args)
      args = set_config args
      config.set_config :single_step, :true
      start_pipe args[0].to_sym
    end

    def start(args)
      args = set_config args
      start_pipe (args[0] || steps.first).to_sym
    end

    def stop(args)
      args = set_config args, :stop
      config.set_config :scheduler, scheduler_type
      cancel_job(args.first == "please")
    end

    def audit(args)
      args = set_config args, :audit
      audit_job(args)
    end

    def clean(args)
      args = set_config args, :clean
      clean_job(args)
    end

    def list_steps(args)
      args = set_config args, :audit
      if args.first == "available"
        puts "Available steps for #{config.script}".red.bold
        self.class.steps.each do |s|
          step = create_step s
          puts "step #{config.step},".green.bold + " runs on #{(step.job_items || [:cohort]).join(".")}".cyan.bold
          puts "  runs tasks #{step.class.tasks.join(", ")}".blue.bold
          puts "  has tasks #{step.class.available_tasks.join(", ")}".yellow.bold if (step.class.available_tasks - step.class.tasks).size > 0
        end
        puts "Available modules for #{config.script}".red.bold
        self.class.modules.each do |m,h|
          puts "module #{m} => ".magenta.bold + h.map{|s,d|
            (self.class.modules[s] ? "#{s}".magenta.bold : "#{s}".green.bold) + (d == true ? "" : "[#{d.join(", ")}]".blue.bold)
          }.join(", ")
        end
      else
        puts "Steps for #{config.script}".red.bold
        steps.each do |s|
          step = create_step s
          puts "step #{config.step},".green.bold + " runs on #{(step.job_items || [:cohort]).join(".")}".cyan.bold
          puts "  runs tasks #{step.tasks.join(", ")}".blue.bold
        end
      end
    end

    def init(args)
      # do something

      if !args.first || !usages[args.first.to_sym]
        usage
        exit
      end

      cmd = args.shift.to_sym

      return if check_usage cmd, args
      send cmd, args
    end
  end
end
