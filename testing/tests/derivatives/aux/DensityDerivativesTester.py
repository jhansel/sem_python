import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from Density import Density, DensityParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

params = DensityParameters()
params.set("p_function", None)
test_aux = Density(params)
test_var = "rho"
params = TestAuxParameters()
params.set("var", "vf")
params.set("other_vars", ["vf1"])
params.set("coefs", [2.0])
vf_aux = TestAux(params)
other_aux = {"vf" : vf_aux}
other_vars = ["vf"]
root_vars = ["vf1", "arho"]

class DensityDerivativesTester(unittest.TestCase):
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
