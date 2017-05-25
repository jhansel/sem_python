import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from SpecificVolume import SpecificVolume, SpecificVolumeParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# specific volume aux
params = SpecificVolumeParameters()
params.set("p_function", None)
test_aux = SpecificVolume(params)
test_var = "v"

# density aux
params = TestAuxParameters()
params.set("var", "rho")
params.set("other_vars", ["vf1", "arho"])
params.set("coefs", [2.0, 3.0])
rho_aux = TestAux(params)

other_aux = {"rho" : rho_aux}
other_vars = ["rho"]
root_vars = ["vf1", "arho"]

class SpecificVolumeDerivativesTester(unittest.TestCase):
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
