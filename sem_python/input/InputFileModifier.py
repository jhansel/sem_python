class InputFileModifier(object):

    def __init__(self):
        self.blocks_to_remove = list()
        self.subblocks_to_remove = list()
        self.block_parameters_to_remove = list()
        self.subblock_parameters_to_remove = list()
        self.block_parameter_changes = list()
        self.subblock_parameter_changes = list()

    def removeBlock(self, block):
        self.blocks_to_remove.append(block)

    def removeSubblock(self, block, subblock):
        self.subblocks_to_remove.append((block, subblock))

    def removeBlockParam(self, block, param):
        self.block_parameters_to_remove.append((block, param))

    def removeSubblockParam(self, block, subblock, param):
        self.subblock_parameters_to_remove.append((block, subblock, param))

    def modifyBlockParam(self, block, param, value):
        self.block_parameter_changes.append((block, param, value))

    def modifySubblockParam(self, block, subblock, param, value):
        self.subblock_parameter_changes.append((block, subblock, param, value))
