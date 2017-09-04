import unittest
import os

import sem
from sem_python.input.ParameterModification import BlockParameterModification
from .. import CSVTester

class SourceSplittingTester(unittest.TestCase):
  def test1Phase(self):
    test_dir = test_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep

    mods = list()
    mods.append(BlockParameterModification("NonlinearSolver", "verbose", False))
    mods.append(BlockParameterModification("Executioner", "verbose", False))
    mods.append(BlockParameterModification("Output", "plot_solution", False))
    mods.append(BlockParameterModification("Output", "save_solution", True))
    mods.append(BlockParameterModification("Output", "solution_file", "solution_1phase.csv"))
    sem.run(test_dir + "source_splitting_1phase.in", mods)

    csv_tester = CSVTester(test_dir, "solution_1phase.csv")
    self.assertTrue(csv_tester.filesAreEqual())
