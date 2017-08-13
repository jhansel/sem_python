import re

from conversion_utilities import stringToInt, stringToFloat
from error_utilities import error, inputError
from file_utilities import checkFileExists

## Class for parsing input files to retrieve parameters
class InputFileParser(object):
  def __init__(self):
    self.level = 0
    self.i_line = 1
    self.block = "None"
    self.subblock = "None"

    self.block_data = dict()
    self.subblock_data = dict()
    self.subblock_list = dict()

    # regular expressions
    self.comment_regex = re.compile(r'\#.*')
    self.block_begin_regex = re.compile(r'\[\w+\]')
    self.block_end_regex = re.compile(r'\[\]')
    self.input_entry_regex = re.compile(r'\w+\s*=\s*\S+')

  ## Parses the input file
  # @param filename  name of the input file to parse
  def parse(self, filename):
    checkFileExists(filename)
    with open(filename) as inp:
      for line in inp:
        self.processLine(line)
        self.i_line += 1
    self.performFinalChecks()

  ## Processes a single line of the input file
  # @param line  line to be processed
  def processLine(self, line):
    # trim off the comment, if any
    line = self.comment_regex.split(line)[0]

    # trim off the leading and trailing whitespace
    line = line.lstrip()
    line = line.rstrip()

    # pass if blank line
    if line == "":
      return
    # check if line is a block begin marker
    elif self.block_begin_regex.match(line):
      self.processBlockBegin(line)
      return
    # check if line is a block end marker
    elif self.block_end_regex.match(line):
      self.processBlockEnd()
      return
    # check if line is a key/value entry
    elif self.input_entry_regex.match(line):
      self.processInputEntry(line)
      return
    else:
      inputError(self.i_line, "Invalid line format.")

  ## Processes a line matching a block begin marker
  # @param line  line to be processed
  def processBlockBegin(self, line):
    # increment block nesting level
    self.level += 1

    # check that sub-block was not nested inside another sub-block
    if self.level > 2:
      inputError(self.i_line, "Blocks cannot be nested within sub-blocks.")

    # get the name of the block or sub-block
    block_or_subblock = line.lstrip("[").rstrip("]")

    # determine whether a block or sub-block was entered
    if self.level == 1:
      entered_block = True
      self.block = block_or_subblock
    else:
      entered_block = False
      self.subblock = block_or_subblock

    # if a block (not a sub-block) was entered
    if entered_block:
      # check that block name is valid
      if self.block in self.block_data:
        inputError(self.i_line, "The block name '" + self.block + "' already exists.")
      else:
        self.block_data[self.block] = dict()
        self.subblock_data[self.block] = dict()
        self.subblock_list[self.block] = list()
    # else a sub-block was entered
    else:
      # check that sub-block name has not already been used
      if self.subblock in self.subblock_data[self.block]:
        inputError(self.i_line, "The sub-block name '" + self.subblock + "' has already been used.")
      else:
        self.subblock_data[self.block][self.subblock] = dict()
        self.subblock_list[self.block].append(self.subblock)

  ## Processes a line matching a block end marker
  def processBlockEnd(self):
    if self.level == 1:
      # reset block name
      self.block = "None"
    elif self.level == 2:
      # reset sub-block name
      self.subblock = "None"
    else:
      inputError(self.i_line, "There is no block or sub-block to end.")

    # update block level
    self.level -= 1

  ## Processes an input entry (key/value pair)
  # @param line  line to be processed
  def processInputEntry(self, line):
    # get key and value
    fields = line.split("=")
    key = fields[0].rstrip()
    value = fields[1].lstrip()

    # store key/value pair in the appropriate dictionary
    if self.level == 1:
      if key in self.block_data[self.block]:
        inputError(self.i_line, "Key '" + key + "' already exists in block.")
      else:
        self.block_data[self.block][key] = value
    elif self.level == 2:
      if key in self.subblock_data[self.block][self.subblock]:
        inputError(self.i_line, "Key '" + key + "' already exists in sub-block.")
      else:
        self.subblock_data[self.block][self.subblock][key] = value
    else:
      inputError(self.i_line, "Input can only be in a block or sub-block.")

  ## Performs final input checks
  def performFinalChecks(self):
    if self.level != 0:
      error("All blocks and sub-blocks must be ended with '[]'.")

  ## Applies a modification to a parameter
  # @param mod  input file modification object
  def applyModification(self, mod):
    block = mod.block
    param = mod.param
    value = mod.value
    self.assertBlockExists(block)
    self.block_data[block][param] = value

  ## Gets dictionary of a block's data
  # @param block  block from which to take data
  def getBlockData(self, block):
    self.assertBlockExists(block)
    return self.block_data[block]

  ## Gets dictionary of a sub-block's data
  # @param block  block to which sub-block belongs
  # @param subblock  sub-block from which to take data
  def getSubblockData(self, block, subblock):
    self.assertSubblockExists(block, subblock)
    return self.subblock_data[block][subblock]

  ## Gets list of names of sub-blocks for a given block in order of occurrence
  # @param block  block for which to get the names of sub-blocks
  def getSubblockNames(self, block):
    self.assertBlockExists(block)
    return self.subblock_list[block]

  ## Checks if a block exists
  # @param[in] block  name of block
  def blockExists(self, block):
    return block in self.block_data

  ## Asserts that a block exists
  # @param block  block for which to assert existence
  def assertBlockExists(self, block):
    if not block in self.block_data:
      error("The block '" + block + "' does not exist.")

  ## Asserts that a sub-block exists
  # @param block  block to which subblock belongs
  # @param subblock  sub-block for which to assert existence
  def assertSubblockExists(self, block, subblock):
    self.assertBlockExists(block)
    if not subblock in self.subblock_data[block]:
      error("The subblock '" + subblock + "' does not exist in block '" + block + "'.")

if __name__ == "__main__":
  parser = InputFileParser()
  parser.parse("../../input/one_phase_ss.in")
  print parser.block_data
