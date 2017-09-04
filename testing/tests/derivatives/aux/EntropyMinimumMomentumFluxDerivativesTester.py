import unittest

from sem_python.aux.EntropyMinimumMomentumFlux import EntropyMinimumMomentumFlux, EntropyMinimumMomentumFluxParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from sem_python.aux.VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumMomentumFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumMomentumFlux(params)

# volume fraction aux
params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_arhou1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5, 1.1, 1.3, 1.6])
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
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2", "grad_vf1", "grad_arho1"])
params.set("coefs", [1.3, 2.6, 3.1, 4.2, 4.0, 2.6, 1.4, 2.1, 1.5])
viscflux_arho_aux = TestAux(params)

other_aux = [vf_aux, visccoef_arhou_aux, rho_aux, grad_u_aux, u_aux, viscflux_arho_aux]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2", "grad_vf1", "grad_arho1", "grad_arhou1"]

class EntropyMinimumMomentumFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
