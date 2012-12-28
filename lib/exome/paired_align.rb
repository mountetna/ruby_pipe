require 'pipeline'
require 'exome/align'
require 'exome/config'

module Exome
  class PairedAlign < Pipeline::Script
    steps :align, :merge, :recal, :hybrid_q_c, :mut_det
  end
end
