import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
test_aux = factory.createObject("BerryInterfacePressureBar", {})

# phase-1 pressure
params = dict()
params["var"] = "p1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.1, 1.2, 1.3, 1.4]
p1_aux = factory.createObject("TestAux", params)

# phase-2 pressure
params = dict()
params["var"] = "p2"
params["other_vars"] = ["aA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.7, 1.9, 2.3, 1.6]
p2_aux = factory.createObject("TestAux", params)

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

other_aux = [p1_aux, p2_aux, z1_aux, z2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class BerryInterfacePressureBarDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 5e-5)
