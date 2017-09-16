class InputFileModifier(object):
  def __init__(self):
    self.subblocks_to_remove = list()
    self.block_parameter_changes = list()
    self.subblock_parameter_changes = list()

  def removeSubblock(self, block, subblock):
    self.subblocks_to_remove((block, subblock))

  def modifyBlockParam(self, block, param, value):
    self.block_parameter_changes.append((block, param, value))

  def modifySubblockParam(self, block, subblock, param, value):
    self.subblock_parameter_changes.append((block, subblock, param, value))
