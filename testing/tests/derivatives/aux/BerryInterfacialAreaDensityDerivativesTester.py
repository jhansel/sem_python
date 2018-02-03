import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# test aux
params = dict()
params["a_int_min"] = 0.1
params["a_int_max"] = 0.7
test_aux = factory.createObject("BerryInterfacialAreaDensity", params)

# phase-1 volume fraction aux
params = dict()
params["var"] = "vf1"
params["other_vars"] = ["aA1"]
params["coefs"] = [1.5]
vf1_aux = factory.createObject("TestAux", params)

other_aux = [vf1_aux]
root_vars = ["aA1"]


class BerryInterfacialAreaDensityDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
