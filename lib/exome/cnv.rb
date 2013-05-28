require 'pipeline'
require 'exome/config'
require 'exome/copy_number'

module Exome
  class Cnv 
    include Pipeline::Script
    runs_steps :copy_number_prep, :copy_number
    #exclude_task :copy_number, :update_mutations
  end
end
