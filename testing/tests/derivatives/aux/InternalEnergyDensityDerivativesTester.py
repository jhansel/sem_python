import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# internal energy density aux
test_aux = factory.createObject("InternalEnergyDensity", {"phase": 0})

# density aux
params = dict()
params["var"] = "rho1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [2.0, 3.0]
rho_aux = factory.createObject("TestAux", params)

# specific internal energy aux
params = dict()
params["var"] = "e1"
params["other_vars"] = ["arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [2.5, 3.5, 4.5]
e_aux = factory.createObject("TestAux", params)

other_aux = [rho_aux, e_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]


class InternalEnergyDensityDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
