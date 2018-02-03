import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester


def computeSoundSpeed(v, e):
    v_slope = 2.0
    e_slope = 3.0
    c = v_slope * v + e_slope * e
    return (c, v_slope, e_slope)

factory = Factory()

# sound speed aux
params = dict()
params["phase"] = 0
params["c_function"] = computeSoundSpeed
test_aux = factory.createObject("SoundSpeed", params)

# specific volume aux
params = dict()
params["var"] = "v1"
params["other_vars"] = ["aA1", "arhoA1"]
params["coefs"] = [2.0, 3.0]
v_aux = factory.createObject("TestAux", params)

# pressure aux
params = dict()
params["var"] = "p1"
params["other_vars"] = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]
params["coefs"] = [1.5, 2.5, 3.5, 4.5]
p_aux = factory.createObject("TestAux", params)

other_aux = [v_aux, p_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]


class SoundSpeedDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
