import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# viscous flux aux
test_aux = factory.createObject("EntropyMinimumMassFlux", {"phase": 0})

# volume fraction aux
vf_aux = factory.createObject("VolumeFractionPhase1", {"phase": 0})

# viscous coefficient aux
params = dict()
params["var"] = "visccoef_arhoA1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.5, 2.5, 3.5, 4.5, 1.1, 1.7, 2.1]
visccoef_arho_aux = factory.createObject("TestAux", params)

# density gradient aux
params = dict()
params["var"] = "grad_rho1"
params["other_vars"] = ["grad_aA1", "grad_arhoA1"]
params["coefs"] = [1.2, 2.4]
grad_rho_aux = factory.createObject("TestAux", params)

# density aux
params = dict()
params["var"] = "rho1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [1.4, 2.3]
rho_aux = factory.createObject("TestAux", params)

# viscous flux aux
params = dict()
params["var"] = "viscflux_aA1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1"]
params["coefs"] = [1.3, 2.6, 3.1, 1.7, 2.1, 1.1, 4.2, 4.0]
viscflux_aA1_aux = factory.createObject("TestAux", params)

other_aux = [vf_aux, visccoef_arho_aux, rho_aux, grad_rho_aux, viscflux_aA1_aux]
root_vars = [
    "aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1", "grad_arhoA1"
]
constant_data = {"A": 0.2}


class EntropyMinimumMassFluxDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data=constant_data)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
