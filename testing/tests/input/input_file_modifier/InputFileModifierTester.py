import unittest

from InputFileModifier import InputFileModifier
from InputFileParser import InputFileParser

class InputFileModifierTester(unittest.TestCase):
  def setUp(self):
    # create modifications
    input_file_modifier = InputFileModifier()
    input_file_modifier.removeBlock("BlockC")
    input_file_modifier.removeSubblock("BlockA", "SubblockAB")
    input_file_modifier.removeBlockParam("BlockA", "paramA3")
    input_file_modifier.removeSubblockParam("BlockA", "SubblockAA", "subparamAA3")
    input_file_modifier.modifyBlockParam("BlockA", "paramA1", "new_value_paramA1")
    input_file_modifier.modifySubblockParam("BlockA", "SubblockAA", "subparamAA1", "new_value_subparamAA1")

    # parse the input file and apply modifications
    input_file = "tests/input/input_file_parser/base.in"
    self.input_file_parser = InputFileParser()
    self.input_file_parser.parse(input_file)
    self.input_file_parser.applyModifications(input_file_modifier)

  def testBlockRemoval(self):
    no_BlockC = ("BlockC" not in self.input_file_parser.block_data) \
      and ("BlockC" not in self.input_file_parser.subblock_data)
    self.assertTrue(no_BlockC)

  def testSubblockRemoval(self):
    no_SubblockAB = ("SubblockAB" not in self.input_file_parser.subblock_data) \
      and ("SubblockAB" not in self.input_file_parser.subblock_list["BlockA"])
    self.assertTrue(no_SubblockAB)

  def testBlockParameterRemoval(self):
    no_paramA3 = "paramA3" not in self.input_file_parser.block_data["BlockA"]
    self.assertTrue(no_paramA3)

  def testSubblockParameterRemoval(self):
    no_subparamAA3 = "subparamAA3" not in self.input_file_parser.subblock_data["BlockA"]["SubblockAA"]
    self.assertTrue(no_subparamAA3)

  def testBlockParameterPreserve(self):
    old_paramA2 = self.input_file_parser.block_data["BlockA"]["paramA2"] == "old_value_paramA2"
    self.assertTrue(old_paramA2)

  def testSubblockParameterPreserve(self):
    old_subparamAA2 = self.input_file_parser.subblock_data["BlockA"]["SubblockAA"]["subparamAA2"] == "old_value_subparamAA2"
    self.assertTrue(old_subparamAA2)

  def testBlockParameterChange(self):
    new_paramA1 = self.input_file_parser.block_data["BlockA"]["paramA1"] == "new_value_paramA1"
    self.assertTrue(new_paramA1)

  def testSubblockParameterChange(self):
    new_subparamAA1 = self.input_file_parser.subblock_data["BlockA"]["SubblockAA"]["subparamAA1"] == "new_value_subparamAA1"
    self.assertTrue(new_subparamAA1)
