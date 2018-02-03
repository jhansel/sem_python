import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
test_aux = factory.createObject("BerryInterfacePressure")

# bar interface pressure aux
params = dict()
params["var"] = "pI_bar"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
pI_bar_aux = factory.createObject("TestAux", params)

# phase-1 velocity aux
params = dict()
params["var"] = "u1"
params["other_vars"] = ["arhoA1", "arhouA1"]
params["coefs"] = [2.3, 3.4]
u1_aux = factory.createObject("TestAux", params)

# phase-2 velocity aux
params = dict()
params["var"] = "u2"
params["other_vars"] = ["arhoA2", "arhouA2"]
params["coefs"] = [2.2, 3.3]
u2_aux = factory.createObject("TestAux", params)

# phase-1 acoustic impedance aux
params = dict()
params["var"] = "z1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.6, 2.3, 4.5, 2.1]
z1_aux = factory.createObject("TestAux", params)

# phase-2 acoustic impedance aux
params = dict()
params["var"] = "z2"
params["other_vars"] = ["aA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.2, 3.2, 4.1, 2.4]
z2_aux = factory.createObject("TestAux", params)

other_aux = [pI_bar_aux, u1_aux, u2_aux, z1_aux, z2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
constant_data_positive = {"grad_aA1": 0.6}
constant_data_negative = {"grad_aA1": -0.6}


class BerryInterfacePressureDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def testPositiveGradient(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data_positive)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)

    def testNegativeGradient(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data_negative)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
