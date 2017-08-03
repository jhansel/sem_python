import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from AcousticImpedance import AcousticImpedance, AcousticImpedanceParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = AcousticImpedanceParameters()
params.set("phase", 0)
test_aux = AcousticImpedance(params)
test_var = "rhoc1"

# density aux
params = TestAuxParameters()
params.set("var", "rho1")
params.set("other_vars", ["vf1", "arho1"])
params.set("coefs", [1.4, 2.5])
rho_aux = TestAux(params)

# sound speed aux
params = TestAuxParameters()
params.set("var", "c1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
c_aux = TestAux(params)

other_aux = {"rho1": rho_aux, "c1": c_aux}
other_vars = ["rho1", "c1"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1"]

class AcousticImpedanceDerivativesTester(unittest.TestCase):
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
