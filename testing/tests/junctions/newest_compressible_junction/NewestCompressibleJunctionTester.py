import unittest

import sem
from InputFileModifier import InputFileModifier
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class NewestCompressibleJunctionTester(unittest.TestCase):
  def runDerivativeTest(self, test_option, use_zero_velocity):
    tester = JunctionTester("NewestCompressibleJunction", rel_tol=5e-5)
    matched = tester.checkJacobian(test_option, use_zero_velocity=use_zero_velocity)
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianWeakNonzeroVelocity(self):
    self.runDerivativeTest("weak", False)

  def testJacobianBothNonzeroVelocity(self):
    self.runDerivativeTest("both", False)

  def testJacobianWeakZeroVelocity(self):
    self.runDerivativeTest("weak", True)

  def testJacobianBothZeroVelocity(self):
    self.runDerivativeTest("both", True)
