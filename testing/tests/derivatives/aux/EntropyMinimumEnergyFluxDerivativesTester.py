import unittest

from EntropyMinimumEnergyFlux import EntropyMinimumEnergyFlux, EntropyMinimumEnergyFluxParameters
from TestAux import TestAux, TestAuxParameters
from VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumEnergyFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumEnergyFlux(params)

# volume fraction aux
params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_arhoE1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5, 1.1, 1.3, 2.1])
visccoef_arhoE_aux = TestAux(params)

# internal energy density gradient aux
params = TestAuxParameters()
params.set("var", "grad_rhoe1")
params.set("other_vars", ["grad_vf1", "grad_arho1", "grad_arhou1", "grad_arhoE1"])
params.set("coefs", [1.3, 2.2, 3.3, 4.0])
grad_rhoe_aux = TestAux(params)

# internal energy density aux
params = TestAuxParameters()
params.set("var", "rhoe1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.4, 2.3, 3.4, 4.1])
rhoe_aux = TestAux(params)

# volume fraction viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_vf1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "grad_vf1"])
params.set("coefs", [1.3, 2.6, 3.1, 4.2, 4.0])
viscflux_vf_aux = TestAux(params)

# velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [1.7, 2.1])
u_aux = TestAux(params)

# mass viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_arho1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2", "grad_vf1", "grad_arho1"])
params.set("coefs", [1.2, 2.5, 3.0, 4.1, 3.9, 2.5, 1.4, 1.3, 1.8])
viscflux_arho_aux = TestAux(params)

# momentum viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_arhou1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2", "grad_vf1", "grad_arho1", "grad_arhou1"])
params.set("coefs", [1.1, 2.2, 3.8, 4.0, 3.7, 2.3, 1.2, 2.1, 1.5, 1.4])
viscflux_arhou_aux = TestAux(params)

other_aux = [vf_aux, visccoef_arhoE_aux, grad_rhoe_aux, rhoe_aux,
             viscflux_vf_aux, u_aux, viscflux_arho_aux, viscflux_arhou_aux]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2", "grad_vf1", "grad_arho1", "grad_arhou1", "grad_arhoE1"]

class EntropyMinimumEnergyFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-5)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
