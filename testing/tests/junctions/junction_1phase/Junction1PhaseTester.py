import unittest

import sem
from InputFileModifier import InputFileModifier
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class Junction1PhaseTester(unittest.TestCase):
  def runDerivativeTest(self, test_option):
    tester = JunctionTester("Junction1Phase")
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
  tester = JunctionTester("Junction1Phase", verbose=True)
  _ = tester.checkJacobian("weak")
  _ = tester.checkJacobian("strong")
