import unittest

import sem
from ParameterModification import BlockParameterModification
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class NewerCompressibleJunctionTester(unittest.TestCase):
  def runDerivativeTest(self, test_option):
    tester = JunctionTester("NewerCompressibleJunction", rel_tol=5e-6)
    matched = tester.checkJacobian(test_option)
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianWeak(self):
    self.runDerivativeTest("weak")

  def testJacobianBoth(self):
    self.runDerivativeTest("both")

if __name__ == "__main__":
  tester = JunctionTester("NewerCompressibleJunction", verbose=True)
  _ = tester.checkJacobian("weak")
  _ = tester.checkJacobian("both")
