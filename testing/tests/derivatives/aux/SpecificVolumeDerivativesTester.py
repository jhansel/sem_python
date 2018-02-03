import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# specific volume aux
test_aux = factory.createObject("SpecificVolume", {"phase": 0})

# density aux
params = dict()
params["var"] = "rho1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [2.0, 3.0]
rho_aux = factory.createObject("TestAux", params)

other_aux = [rho_aux]
root_vars = ["aA1", "arhoA1"]


class SpecificVolumeDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
