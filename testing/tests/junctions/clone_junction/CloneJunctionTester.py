import unittest

from InputFileModifier import InputFileModifier
from SolutionTester import SolutionTester
from JunctionTester import JunctionTester

class CloneJunctionTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "tests/junctions/clone_junction/"
    input_file = "../problems/junction/junction.in"

    input_file_modifier = InputFileModifier()
    input_file_modifier.modifySubblockParam("Junctions", "junction1", "type", "CloneJunction")
    input_file_modifier.removeSubblockParam("Junctions", "junction1", "phase")
    input_file_modifier.modifyBlockParam("Executioner", "end_time", 0.05)

    solution_tester = SolutionTester(test_dir, input_file, input_file_modifier)
    self.assertTrue(solution_tester.solutionsAreEqual())

  def runDerivativeTest(self, test_option):
    tester = JunctionTester("CloneJunction")
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
  tester = JunctionTester("CloneJunction", verbose=True)
  _ = tester.checkJacobian("weak")
  _ = tester.checkJacobian("strong")
