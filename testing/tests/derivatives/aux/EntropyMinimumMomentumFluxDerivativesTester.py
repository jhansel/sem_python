import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# viscous flux aux
test_aux = factory.createObject("EntropyMinimumMomentumFlux", {"phase": 0})

# volume fraction aux
vf_aux = factory.createObject("VolumeFractionPhase1", {"phase": 0})

# viscous coefficient aux
params = dict()
params["var"] = "visccoef_arhouA1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.5, 2.5, 3.5, 4.5, 1.1, 1.3, 1.6]
visccoef_arhouA_aux = factory.createObject("TestAux", params)

# density aux
params = dict()
params["var"] = "rho1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [1.4, 2.3]
rho_aux = factory.createObject("TestAux", params)

# velocity gradient aux
params = dict()
params["var"] = "grad_u1"
params["other_vars"] = ["grad_arhoA1", "grad_arhouA1"]
params["coefs"] = [1.2, 2.4]
grad_u_aux = factory.createObject("TestAux", params)

# velocity aux
params = dict()
params["var"] = "u1"
params["other_vars"] = ["arhoA1", "arhouA1"]
params["coefs"] = [1.7, 2.1]
u_aux = factory.createObject("TestAux", params)

# viscous flux aux
params = dict()
params["var"] = "viscflux_arhoA1"
params[
    "other_vars"] = [
        "aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1",
        "grad_arhoA1"
    ]
params["coefs"] = [1.3, 2.6, 3.1, 4.2, 4.0, 2.6, 1.4, 2.1, 1.5]
viscflux_arho_aux = factory.createObject("TestAux", params)

other_aux = [vf_aux, visccoef_arhouA_aux, rho_aux, grad_u_aux, u_aux, viscflux_arho_aux]
root_vars = [
    "aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1",
    "grad_arhoA1", "grad_arhouA1"
]
constant_data = {"A": 0.6}


class EntropyMinimumMomentumFluxDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data=constant_data)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
