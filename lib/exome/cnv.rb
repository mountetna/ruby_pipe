require 'pipeline'
require 'exome/config'
require 'exome/copy_number'

module Exome
  class Cnv 
    include Pipeline::Script
    runs_steps :copy_number_prep, :copy_number
    def_module :default, { :copy_number_prep => true, :copy_number => true }
  end
end
