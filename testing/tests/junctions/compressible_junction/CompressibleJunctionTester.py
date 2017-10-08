import unittest

import sem
from sem_python.input.InputFileModifier import InputFileModifier
from ....src.testers.CSVTester import CSVTester
from ....src.testers.JunctionTester import JunctionTester

class CompressibleJunctionTester(unittest.TestCase):
  def runDerivativeTest(self, test_option, use_momentum_flux_balance, use_lm):
    tester = JunctionTester("CompressibleJunction")
    matched = tester.checkJacobian(test_option,
      junction_params={"use_momentum_flux_balance": use_momentum_flux_balance, "use_lm": use_lm})
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianTechnique1NoLM(self):
    self.runDerivativeTest("both", False, False)

  def testJacobianTechnique2NoLM(self):
    self.runDerivativeTest("both", True, False)

if __name__ == "__main__":
  tester = JunctionTester("CompressibleJunction", verbose=True)
  print "\nDirect, Stagnation pressure"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": False, "use_lm": False})
  print "\nDirect, Momentum flux balance"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": True, "use_lm": False})
  print "\nLagrange Multiplier, Stagnation pressure"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": False, "use_lm": True})
  print "\nLagrange Multiplier, Momentum flux balance"
  _ = tester.checkJacobian("both", junction_params={"use_momentum_flux_balance": True, "use_lm": True})
