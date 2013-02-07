require 'fileutils'
module Pipeline
  module Script
    include Pipeline::Base
    include Pipeline::Logger

    module ClassMethods 
      attr_reader :steps
      def runs_steps(*step_list)
        @steps = step_list
      end
    end

    def self.included(base)
      base.extend(ClassMethods) if ClassMethods
    end

    class_var :steps

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
      self.class.sister_class(s).new(self)
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

    def usage(cmd=nil)
      puts "Possible commands:"
      cmds = {
        :start => {
          :cmd => "start <config_file.yml> [<step>]" ,
          :blurb => "Start the pipeline at the beginning or at <step>"
        },
        :run_step => {
          :cmd => "run_step <config_file.yml> <step_name>",
          :blurb => "Run just the named step"
        },
        :audit => {
          :cmd => "audit <config_file.yml>",
          :blurb => "Audit the pipeline to see which steps are complete."
        },
        :generate => {
          :cmd => "generate <job_name>",
          :blurb => "Generate a new config file."
        }
      }
      if cmds[cmd]
        puts " %-50s" % "#{cmds[cmd][:cmd]}" + "# #{cmds[cmd][:blurb]}".cyan
      else
        cmds.each do |a,c|
          puts " %-50s" % "#{c[:cmd]}" + "# #{c[:blurb]}".cyan
        end
      end
      exit
    end

    def generate(args)
      if args.length < 1
        usage :generate
      end
      self.class.daughter_class(:config_generator).new *args
    end

    def init(args)
      # do something
      cmd = args.shift

      ENV['CONFIG'] = args.shift if [ "start", "run_step", "audit" ].include? cmd
      
      case cmd
      when "start"
        config.set_config :action, :init
        start_pipe (args[0] || steps.first).to_sym
      when "run_step"
        exit unless args.length == 1
        config.set_config :action, :init
        config.set_config :single_step, :true
        start_pipe args[0].to_sym
      when "audit"
        config.set_config :action, :init
        audit(args)
      when "generate"
        generate(args)
      else
        usage
      end
    end
  end
end
