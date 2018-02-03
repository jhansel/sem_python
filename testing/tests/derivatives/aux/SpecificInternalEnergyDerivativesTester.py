import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# specific internal energy aux
test_aux = factory.createObject("SpecificInternalEnergy", {"phase": 0})

# velocity aux
params = dict()
params["var"] = "u1"
params["other_vars"] = ["arhoA1", "arhouA1"]
params["coefs"] = [1.1, 1.2]
u_aux = factory.createObject("TestAux", params)

# specific total energy aux
params = dict()
params["var"] = "E1"
params["other_vars"] = ["arhoA1", "arhoEA1"]
params["coefs"] = [3.3, 3.9]
E_aux = factory.createObject("TestAux", params)

other_aux = [u_aux, E_aux]
root_vars = ["arhoA1", "arhouA1", "arhoEA1"]


class SpecificInternalEnergyDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
