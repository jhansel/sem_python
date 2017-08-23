import unittest

import sem
from InputFileModification import InputFileModification
from CSVTester import CSVTester

class TwoPhaseNoninteractingSteadyStateTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/solution/two_phase_noninteracting_ss/"
    sem.run(test_dir + "two_phase_noninteracting_ss.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(InputFileModification("NonlinearSolver", "verbose", True))
  mods.append(InputFileModification("Output", "save_solution", False))
  sem.run("two_phase_noninteracting_ss.in", mods)
