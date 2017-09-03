import unittest

from EntropyMinimumMassFlux import EntropyMinimumMassFlux, EntropyMinimumMassFluxParameters
from TestAux import TestAux, TestAuxParameters
from VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumMassFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumMassFlux(params)

# volume fraction aux
params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_arhoA1")
params.set("other_vars", ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5, 1.1, 1.7, 2.1])
visccoef_arho_aux = TestAux(params)

# density gradient aux
params = TestAuxParameters()
params.set("var", "grad_rho1")
params.set("other_vars", ["grad_vf1", "grad_arhoA1"])
params.set("coefs", [1.2, 2.4])
grad_rho_aux = TestAux(params)

# density aux
params = TestAuxParameters()
params.set("var", "rho1")
params.set("other_vars", ["vf1", "arhoA1"])
params.set("coefs", [1.4, 2.3])
rho_aux = TestAux(params)

# viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_vf1")
params.set("other_vars", ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_vf1"])
params.set("coefs", [1.3, 2.6, 3.1, 1.7, 2.1, 1.1, 4.2, 4.0])
viscflux_vf1_aux = TestAux(params)

other_aux = [vf_aux, visccoef_arho_aux, rho_aux, grad_rho_aux, viscflux_vf1_aux]
root_vars = ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_vf1", "grad_arhoA1"]

class EntropyMinimumMassFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
