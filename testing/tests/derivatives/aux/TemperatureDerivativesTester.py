import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester


def computeTemperature(v, e):
    v_slope = 2.0
    e_slope = 3.0
    T = v_slope * v + e_slope * e
    return (T, v_slope, e_slope)

factory = Factory()

# temperature aux
params = dict()
params["phase"] = 0
params["T_function"] = computeTemperature
test_aux = factory.createObject("Temperature", params)

# specific volume aux
params = dict()
params["var"] = "v1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [2.0, 3.0]
v_aux = factory.createObject("TestAux", params)

# specific internal energy aux
params = dict()
params["var"] = "e1"
params["other_vars"] = ["arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [2.5, 3.5, 4.5]
e_aux = factory.createObject("TestAux", params)

other_aux = [v_aux, e_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]


class TemperatureDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
