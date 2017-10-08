import unittest

from InputFileModifier import InputFileModifier
from ....src.testers.SolutionTester import SolutionTester
from ....src.testers.JunctionTester import JunctionTester

class EqualSolutionLM1PhaseJunctionTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "testing/tests/junctions/equal_solution_lm_1phase_junction/"
    input_file = "problems/junction/junction.in"

    input_file_modifier = InputFileModifier()
    input_file_modifier.modifySubblockParam("Junctions", "junction1", "type", "EqualSolutionLM1PhaseJunction")
    input_file_modifier.modifyBlockParam("Executioner", "end_time", 0.05)

    solution_tester = SolutionTester(test_dir, input_file, input_file_modifier)
    self.assertTrue(solution_tester.solutionsAreEqual())

  def runDerivativeTest(self, test_option):
    tester = JunctionTester("EqualSolutionLM1PhaseJunction")
    matched = tester.checkJacobian(test_option)
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianWeak(self):
    self.runDerivativeTest("weak")

  def testJacobianStrong(self):
    self.runDerivativeTest("strong")

if __name__ == "__main__":
  tester = JunctionTester("EqualSolutionLM1PhaseJunction", verbose=True)
  _ = tester.checkJacobian("weak")
  _ = tester.checkJacobian("strong")
