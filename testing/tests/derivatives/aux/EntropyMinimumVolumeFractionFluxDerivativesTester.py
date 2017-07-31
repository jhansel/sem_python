import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from EntropyMinimumVolumeFractionFlux import EntropyMinimumVolumeFractionFlux, EntropyMinimumVolumeFractionFluxParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumVolumeFractionFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumVolumeFractionFlux(params)
test_var = "viscflux_vf1"

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_vf1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
visccoef_vf_aux = TestAux(params)

# volume fraction gradient aux
params = TestAuxParameters()
params.set("var", "grad_vf")
params.set("other_vars", ["grad_vf1"])
params.set("coefs", [1.0])
grad_vf_aux = TestAux(params)

other_aux = {"visccoef_vf1" : visccoef_vf_aux, "grad_vf" : grad_vf_aux}
other_vars = ["visccoef_vf1", "grad_vf"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "grad_vf1"]

class EntropyMinimumVolumeFractionFluxDerivativesTester(unittest.TestCase):
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
