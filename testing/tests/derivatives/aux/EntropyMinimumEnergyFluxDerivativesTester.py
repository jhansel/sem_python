import unittest

from sem_python.aux.EntropyMinimumEnergyFlux import EntropyMinimumEnergyFlux, EntropyMinimumEnergyFluxParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from sem_python.aux.VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

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
params.set("var", "visccoef_arhoEA1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5, 1.1, 1.3, 2.1])
visccoef_arhoEA_aux = TestAux(params)

# internal energy density gradient aux
params = TestAuxParameters()
params.set("var", "grad_rhoe1")
params.set("other_vars", ["grad_aA1", "grad_arhoA1", "grad_arhouA1", "grad_arhoEA1"])
params.set("coefs", [1.3, 2.2, 3.3, 4.0])
grad_rhoe_aux = TestAux(params)

# internal energy density aux
params = TestAuxParameters()
params.set("var", "rhoe1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.4, 2.3, 3.4, 4.1])
rhoe_aux = TestAux(params)

# volume fraction viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_aA1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "grad_aA1"])
params.set("coefs", [1.3, 2.6, 3.1, 4.2, 4.0])
viscflux_aA_aux = TestAux(params)

# velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arhoA1", "arhouA1"])
params.set("coefs", [1.7, 2.1])
u_aux = TestAux(params)

# mass viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_arhoA1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1", "grad_arhoA1"])
params.set("coefs", [1.2, 2.5, 3.0, 4.1, 3.9, 2.5, 1.4, 1.3, 1.8])
viscflux_arho_aux = TestAux(params)

# momentum viscous flux aux
params = TestAuxParameters()
params.set("var", "viscflux_arhouA1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1", "grad_arhoA1", "grad_arhouA1"])
params.set("coefs", [1.1, 2.2, 3.8, 4.0, 3.7, 2.3, 1.2, 2.1, 1.5, 1.4])
viscflux_arhouA_aux = TestAux(params)

other_aux = [vf_aux, visccoef_arhoEA_aux, grad_rhoe_aux, rhoe_aux,
             viscflux_aA_aux, u_aux, viscflux_arho_aux, viscflux_arhouA_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1", "grad_arhoA1", "grad_arhouA1", "grad_arhoEA1"]

class EntropyMinimumEnergyFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data={"A": 0.3})
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-5)
