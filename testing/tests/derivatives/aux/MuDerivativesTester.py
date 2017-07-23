import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from Mu import Mu, MuParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

def computeMu(T1, T2, beta):
  T1_slope = 2.0
  T2_slope = 3.0
  beta_slope = 4.0
  mu = T1_slope * T1 + T2_slope * T2 + beta_slope * beta
  return (mu, T1_slope, T2_slope, beta_slope)

# interface pressure aux
params = MuParameters()
params.set("mu_function", computeMu)
test_aux = Mu(params)
test_var = "mu"

# phase-1 temperature aux
params = TestAuxParameters()
params.set("var", "T1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
T1_aux = TestAux(params)

# phase-2 temperature aux
params = TestAuxParameters()
params.set("var", "T2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
T2_aux = TestAux(params)

# beta aux
params = TestAuxParameters()
params.set("var", "beta")
params.set("other_vars", ["arho1", "arho2"])
params.set("coefs", [1.7, 2.7])
beta_aux = TestAux(params)

other_aux = {"T1" : T1_aux, "T2" : T2_aux, "beta" : beta_aux}
other_vars = ["T1", "T2", "beta"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class MuDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux, other_vars, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux, other_vars, root_vars)
