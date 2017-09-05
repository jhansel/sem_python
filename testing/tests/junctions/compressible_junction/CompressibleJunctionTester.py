import unittest

import sem
from ParameterModification import BlockParameterModification
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class CompressibleJunctionTester(unittest.TestCase):
  def runDerivativeTest(self, test_option, use_momentum_flux_balance, use_lm):
    tester = JunctionTester("CompressibleJunction")
    matched = tester.checkJacobian(test_option,
      junction_params={"use_momentum_flux_balance": use_momentum_flux_balance, "use_lm": use_lm})
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianTechnique1(self):
    self.runDerivativeTest("both", False)

  def testJacobianTechnique2(self):
    self.runDerivativeTest("both", True)

if __name__ == "__main__":
  tester = JunctionTester("CompressibleJunction", verbose=True)
  print "\nDirect, Stagnation pressure"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": False, "use_lm": False})
  print "\nDirect, Momentum flux balance"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": True, "use_lm": False})
  print "\nLagrange Multiplier, Momentum flux balance"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": True, "use_lm": True})
