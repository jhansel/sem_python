import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
test_aux = factory.createObject("AcousticImpedance", {"phase": 0})

# density aux
params = dict()
params["var"] = "rho1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [1.4, 2.5]
rho_aux = factory.createObject("TestAux", params)

# sound speed aux
params = dict()
params["var"] = "c1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.2, 2.2, 3.2, 4.2]
c_aux = factory.createObject("TestAux", params)

other_aux = [rho_aux, c_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]


class AcousticImpedanceDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
