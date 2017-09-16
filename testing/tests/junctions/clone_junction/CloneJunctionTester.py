import unittest

from ParameterModification import BlockParameterModification, SubblockParameterModification
from SolutionTester import SolutionTester
from JunctionTester import JunctionTester

class CloneJunctionTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "tests/junctions/clone_junction/"
    input_file = "tests/junctions/junction.in"

    mods = list()
    mods.append(SubblockParameterModification("Junctions", "junction1", "type", "CloneJunction"))
    mods.append(BlockParameterModification("Executioner", "end_time", 0.05))

    solution_tester = SolutionTester(test_dir, input_file, mods)
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
