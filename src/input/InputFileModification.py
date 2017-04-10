## Class for modifying an input file value
class InputFileModification(object):
  ## Constructor
  # @param block  block in which parameter to modify resides
  # @param param  parameter to modify
  # @param value  new value of parameter
  def __init__(self, block, param, value):
    self.block = block
    self.param = param
    self.value = value
