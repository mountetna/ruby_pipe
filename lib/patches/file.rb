#!/usr/bin/env ruby

class File
  def self.human_size file
    return 0 if !exists?(file) || !size?(file)
    sz = size file
    case sz
    when 0..1023 
      "%d b" % sz
    when 1024..1024*1024-1 
      "%.1f Kb" % (sz / 1024)
    when 1024*1024..1024*1024*1024-1 
      "%.1f Mb" % (sz / 1024 / 1024)
    when 1024*1024*1024..Inf 
      "%.1f Gb" % (sz / 1024 / 1024 / 1024)
    end
  end
end
