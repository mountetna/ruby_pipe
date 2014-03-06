module Pipeline
  module Usage
    module ClassMethods
      attr_reader :usages

      def usage cmd, expln
        @usages ||= { }
        cmd,args = cmd.split(/ /,2)
        @usages.update Hash[ cmd.to_sym, [ args, required_args(args), expln ] ]
      end

      def required_args args
        args ? args.gsub(/\[.*?\]/,"").scan(/<.*?>/).size : 0
      end
    end

    def self.included(base)
      base.extend ClassMethods
    end

    def usages
      self.class.ancestors.map{ |c|
        c.respond_to?(:usages) ? c.usages : nil
      }.compact.reverse.reduce(:merge) || {}
    end

    def print_usage c, u
      cmd = [ c, u.first ].compact.join " "
      puts " %-50s" % cmd.bold + "# #{u.last}".cyan
    end

    def check_usage cmd, args
      args ||= []
      required = usages[cmd][1]
      if args.size < required
        usage cmd
        return true
      end
    end

    def usage(cmd=nil)
      puts "Commands:"
      if cmd && usages[cmd]
        print_usage cmd, usages[cmd]
      else
        usages.each do |c,u|
          print_usage c, u
        end
      end
      nil
    end
  end
end
