import unittest
import os

import sem
from sem_python.input.ParameterModification import BlockParameterModification
from .. import CSVTester

class LaxFriedrichs2PhaseTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
    sem.run(test_dir + "lax_friedrichs_2phase.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(BlockParameterModification("NonlinearSolver", "verbose", True))
  mods.append(BlockParameterModification("Executioner", "verbose", True))
  mods.append(BlockParameterModification("Output", "save_solution", True))
  mods.append(BlockParameterModification("Output", "solution_file", "solution.csv"))
  sem.run("lax_friedrichs_2phase.in", mods)
