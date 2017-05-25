import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from SpecificInternalEnergy import SpecificInternalEnergy, SpecificInternalEnergyParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# specific internal energy aux
params = SpecificInternalEnergyParameters()
params.set("p_function", None)
test_aux = SpecificInternalEnergy(params)
test_var = "e"

# velocity aux
params = TestAuxParameters()
params.set("var", "u")
params.set("other_vars", ["arho", "arhou"])
params.set("coefs", [2.0, 3.0])
u_aux = TestAux(params)

# specific total energy aux
params = TestAuxParameters()
params.set("var", "E")
params.set("other_vars", ["arho", "arhoE"])
params.set("coefs", [2.5, 3.5])
E_aux = TestAux(params)

other_aux = {"u" : u_aux, "E" : E_aux}
other_vars = ["u", "E"]
root_vars = ["arho", "arhou", "arhoE"]

class SpecificInternalEnergyDerivativesTester(unittest.TestCase):
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
