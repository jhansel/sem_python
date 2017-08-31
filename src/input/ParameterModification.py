## Class for modifying a block parameter
class BlockParameterModification(object):
  ## Constructor
  # @param block  block in which parameter to modify resides
  # @param param  parameter to modify
  # @param value  new value of parameter
  def __init__(self, block, param, value):
    self.is_subblock_param = False
    self.block = block
    self.param = param
    self.value = value

## Class for modifying a sub-block parameter
class SubblockParameterModification(BlockParameterModification):
  ## Constructor
  # @param block  block in which parameter to modify resides
  # @param subblock  sub-block in which parameter to modify resides
  # @param param  parameter to modify
  # @param value  new value of parameter
  def __init__(self, block, subblock, param, value):
    BlockParameterModification.__init__(self, block, param, value)
    self.subblock = subblock
    self.is_subblock_param = True
