import unittest
import os

import sem
from sem_python.input.ParameterModification import BlockParameterModification
from .. import CSVTester

class TwoPhaseNoninteractingSteadyStateTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
    sem.run(test_dir + "two_phase_noninteracting_ss.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(BlockParameterModification("NonlinearSolver", "verbose", True))
  mods.append(BlockParameterModification("Output", "save_solution", False))
  sem.run("two_phase_noninteracting_ss.in", mods)
