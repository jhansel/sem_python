import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
test_aux = factory.createObject("BerryInterfaceVelocityBar", {})

# phase-1 velocity
params = dict()
params["var"] = "u1"
params["other_vars"] = ["arhoA1", "arhouA1"]
params["coefs"] = [1.1, 1.2]
u1_aux = factory.createObject("TestAux", params)

# phase-2 velocity
params = dict()
params["var"] = "u2"
params["other_vars"] = ["arhoA2", "arhouA2"]
params["coefs"] = [1.7, 1.9]
u2_aux = factory.createObject("TestAux", params)

# phase-1 acoustic impedance aux
params = dict()
params["var"] = "z1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [2.1, 2.3, 3.4, 4.5]
z1_aux = factory.createObject("TestAux", params)

# phase-2 acoustic impedance aux
params = dict()
params["var"] = "z2"
params["other_vars"] = ["aA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [2.4, 2.2, 3.3, 4.4]
z2_aux = factory.createObject("TestAux", params)

other_aux = [u1_aux, u2_aux, z1_aux, z2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class BerryInterfaceVelocityBarDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-5)
