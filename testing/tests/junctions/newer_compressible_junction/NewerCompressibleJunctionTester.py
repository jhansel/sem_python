import unittest

import sem
from sem_python.input.InputFileModifier import InputFileModifier
from ....src.testers.CSVTester import CSVTester
from ....src.testers.JunctionTester import JunctionTester

class NewerCompressibleJunctionTester(unittest.TestCase):
  def runDerivativeTest(self, test_option):
    tester = JunctionTester("NewerCompressibleJunction", rel_tol=1e-5)
    matched = tester.checkJacobian(test_option, junction_params={"loss_coefficients": [0.3, 0.4]})
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianWeak(self):
    self.runDerivativeTest("weak")

  def testJacobianBoth(self):
    self.runDerivativeTest("both")
