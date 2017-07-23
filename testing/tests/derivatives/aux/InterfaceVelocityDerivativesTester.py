import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from InterfaceVelocity import InterfaceVelocity, InterfaceVelocityParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

def computeInterfaceVelocity(u1, u2, beta):
  u1_slope = 2.0
  u2_slope = 3.0
  beta_slope = 4.0
  uI = u1_slope * u1 + u2_slope * u2 + beta_slope * beta
  return (uI, u1_slope, u2_slope, beta_slope)

# interface velocity aux
params = InterfaceVelocityParameters()
params.set("uI_function", computeInterfaceVelocity)
test_aux = InterfaceVelocity(params)
test_var = "uI"

# phase-1 velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [1.2, 2.2])
u1_aux = TestAux(params)

# phase-2 velocity aux
params = TestAuxParameters()
params.set("var", "u2")
params.set("other_vars", ["arho2", "arhou2"])
params.set("coefs", [1.5, 2.5])
u2_aux = TestAux(params)

# beta aux
params = TestAuxParameters()
params.set("var", "beta")
params.set("other_vars", ["arho1", "arho2"])
params.set("coefs", [1.7, 3.2])
beta_aux = TestAux(params)

other_aux = {"u1" : u1_aux, "u2" : u2_aux, "beta" : beta_aux}
other_vars = ["u1", "u2", "beta"]
root_vars = ["arho1", "arhou1", "arho2", "arhou2"]

class InterfaceVelocityDerivativesTester(unittest.TestCase):
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
