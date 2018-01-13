import unittest

from sem_python.aux.Pressure import Pressure, PressureParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester


def computePressure(v, e):
    v_slope = 2.0
    e_slope = 3.0
    p = v_slope * v + e_slope * e
    return (p, v_slope, e_slope)


# pressure aux
params = PressureParameters()
params.set("phase", 0)
params.set("p_function", computePressure)
test_aux = Pressure(params)

# specific volume aux
params = TestAuxParameters()
params.set("var", "v1")
params.set("other_vars", ["aA1", "arhoA1"])
params.set("coefs", [2.0, 3.0])
v_aux = TestAux(params)

# specific internal energy aux
params = TestAuxParameters()
params.set("var", "e1")
params.set("other_vars", ["arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [2.5, 3.5, 4.5])
e_aux = TestAux(params)

other_aux = [v_aux, e_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]


class PressureDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
