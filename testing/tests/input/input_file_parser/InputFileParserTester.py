import unittest

from sem_python.input.InputFileParser import InputFileParser

class InputFileParserTester(unittest.TestCase):
  def setUp(self):
    differential_input_file = "testing/tests/input/input_file_parser/differential.in"
    input_file_parser_differential = InputFileParser()
    input_file_parser_differential.parse(differential_input_file)

    # get the name of the base input file and parse it
    base_input_file = input_file_parser_differential.getBlockData("BaseInputFile")["base"]
    self.input_file_parser = InputFileParser()
    self.input_file_parser.parse(base_input_file)

    # apply modifications to base input file parser
    self.input_file_parser.applyDifferentialInputFileParser(input_file_parser_differential)

  def testBlockRemoval(self):
    no_BlockC = ("BlockC" not in self.input_file_parser.block_data) \
      and ("BlockC" not in self.input_file_parser.subblock_data)
    self.assertTrue(no_BlockC)

  def testBlockReplacement(self):
    no_SubblockBA = ("SubblockBA" not in self.input_file_parser.subblock_data["BlockB"]) \
      and ("SubblockBA" not in self.input_file_parser.subblock_list["BlockB"])
    no_paramB1 = "paramB1" not in self.input_file_parser.block_data["BlockB"]
    new_paramB2 = self.input_file_parser.block_data["BlockB"]["paramB2"] == "new_value_paramB2"
    new_subparamBB1 = self.input_file_parser.subblock_data["BlockB"]["SubblockBB"]["subparamBB1"] == "new_value_subparamBB1"
    self.assertTrue(no_SubblockBA and no_paramB1 and new_paramB2 and new_subparamBB1)

  def testSubblockRemoval(self):
    no_SubblockAB = ("SubblockAB" not in self.input_file_parser.subblock_data) \
      and ("SubblockAB" not in self.input_file_parser.subblock_list["BlockA"])
    self.assertTrue(no_SubblockAB)

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
