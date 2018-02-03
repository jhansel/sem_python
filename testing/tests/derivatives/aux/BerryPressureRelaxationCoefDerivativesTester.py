import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
test_aux = factory.createObject("BerryPressureRelaxationCoef")

# interfacial area density aux
params = dict()
params["var"] = "a_int"
params["other_vars"] = ["aA1"]
params["coefs"] = [1.5]
a_int_aux = factory.createObject("TestAux", params)

# phase-1 acoustic impedance aux
params = dict()
params["var"] = "z1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.2, 2.3, 3.4, 4.5]
z1_aux = factory.createObject("TestAux", params)

# phase-2 acoustic impedance aux
params = dict()
params["var"] = "z2"
params["other_vars"] = ["aA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.1, 2.2, 3.3, 4.4]
z2_aux = factory.createObject("TestAux", params)

other_aux = [a_int_aux, z1_aux, z2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class BerryPressureRelaxationCoefDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-6)
