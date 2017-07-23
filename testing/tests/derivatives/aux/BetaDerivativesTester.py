import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from Beta import Beta, BetaParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

def computeBeta(arho1, arho2):
  arho1_slope = 2.0
  arho2_slope = 3.0
  beta = arho1_slope * arho1 + arho2_slope * arho2
  return (beta, arho1_slope, arho2_slope)

# beta aux
params = BetaParameters()
params.set("beta_function", computeBeta)
test_aux = Beta(params)
test_var = "beta"

other_aux = {}
other_vars = []
root_vars = ["arho1", "arho2"]

class BetaDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux, other_vars, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux, other_vars, root_vars)
