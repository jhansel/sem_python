import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# Lax-Friedrichs coefficient aux
test_aux = factory.createObject("LaxFriedrichsCoefficientVolumeFraction", {})

# constant data
constant_data = {"dx": 0.5}

# positive interfacial velocity aux
params = dict()
params["var"] = "uI"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
uI_positive_aux = factory.createObject("TestAux", params)

# negative interfacial velocity aux
params = dict()
params["var"] = "uI"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]
params["coefs"] = [-1.3, -1.5, -1.7, -1.9, -2.1, -2.3, -2.5]
uI_negative_aux = factory.createObject("TestAux", params)

other_aux_positive = [uI_positive_aux]
other_aux_negative = [uI_negative_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]


class LaxFriedrichsCoefficientVolumeFractionDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def testPositiveVelocity(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux_positive, root_vars, constant_data)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)

    def testNegativeVelocity(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux_negative, root_vars, constant_data)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
