import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from InterfacePressure import InterfacePressure, InterfacePressureParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

def computeInterfacePressure(p1, p2, mu):
  p1_slope = 2.0
  p2_slope = 3.0
  mu_slope = 4.0
  pI = p1_slope * p1 + p2_slope * p2 + mu_slope * mu
  return (pI, p1_slope, p2_slope, mu_slope)

# interface pressure aux
params = InterfacePressureParameters()
params.set("pI_function", computeInterfacePressure)
test_aux = InterfacePressure(params)
test_var = "pI"

# phase-1 pressure aux
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
p1_aux = TestAux(params)

# phase-2 pressure aux
params = TestAuxParameters()
params.set("var", "p2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
p2_aux = TestAux(params)

# mu aux
params = TestAuxParameters()
params.set("var", "mu")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.2, 1.5, 1.7, 2.2, 2.5, 2.7, 3.2])
mu_aux = TestAux(params)

other_aux = {"p1" : p1_aux, "p2" : p2_aux, "mu" : mu_aux}
other_vars = ["p1", "p2", "mu"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class InterfacePressureDerivativesTester(unittest.TestCase):
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
