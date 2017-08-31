import unittest

import sem
from ParameterModification import BlockParameterModification
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class FreeBCJunctionTester(unittest.TestCase):
  def runDerivativeTest(self, test_option):
    tester = JunctionTester("FreeBCJunction")
    rel_diffs = tester.checkJacobian(test_option)
    n_i, n_j = rel_diffs.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)

  def testJacobianWeak(self):
    self.runDerivativeTest("weak")

  def testJacobianStrong(self):
    self.runDerivativeTest("strong")

if __name__ == "__main__":
  tester = JunctionTester("FreeBCJunction", verbose=True)
  _ = tester.checkJacobian("weak")
  _ = tester.checkJacobian("strong")
