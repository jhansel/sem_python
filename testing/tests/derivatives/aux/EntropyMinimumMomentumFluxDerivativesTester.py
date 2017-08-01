import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from EntropyMinimumMomentumFlux import EntropyMinimumMomentumFlux, EntropyMinimumMomentumFluxParameters
from TestAux import TestAux, TestAuxParameters
from VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumMomentumFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumMomentumFlux(params)
test_var = "viscflux_arhou1"

# volume fraction aux
params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_arhou1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
visccoef_arhou_aux = TestAux(params)

# density aux
params = TestAuxParameters()
params.set("var", "rho1")
params.set("other_vars", ["vf1", "arho1"])
params.set("coefs", [1.4, 2.3])
rho_aux = TestAux(params)

# velocity gradient aux
params = TestAuxParameters()
params.set("var", "grad_u1")
params.set("other_vars", ["grad_arho1", "grad_arhou1"])
params.set("coefs", [1.2, 2.4])
grad_u_aux = TestAux(params)

# velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [1.7, 2.1])
u_aux = TestAux(params)

# viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_arho1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "grad_vf1", "grad_arho1"])
params.set("coefs", [1.3, 2.6, 3.1, 4.2, 4.0, 2.6])
viscflux_arho_aux = TestAux(params)

other_aux = {"vf1" : vf_aux, "visccoef_arhou1" : visccoef_arhou_aux, "rho1" : rho_aux,
  "grad_u1" : grad_u_aux, "u1" : u_aux, "viscflux_arho1" : viscflux_arho_aux}
other_vars = ["vf1", "visccoef_arhou1", "rho1", "grad_u1", "u1", "viscflux_arho1"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "grad_vf1", "grad_arho1", "grad_arhou1"]

class EntropyMinimumMomentumFluxDerivativesTester(unittest.TestCase):
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
