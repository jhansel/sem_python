import unittest

import sem
from ParameterModification import BlockParameterModification
from CSVTester import CSVTester

class OnePhaseSteadyStateTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/solution/one_phase_ss/"
    sem.run(test_dir + "one_phase_ss.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(BlockParameterModification("NonlinearSolver", "verbose", True))
  mods.append(BlockParameterModification("Output", "save_solution", False))
  sem.run("one_phase_ss.in", mods)
