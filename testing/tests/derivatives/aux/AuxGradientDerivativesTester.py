import unittest

from sem_python.aux.AuxGradient import AuxGradient, AuxGradientParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# gradient of test aux
params = AuxGradientParameters()
params.set("aux", "testaux")
params.set("variable_names", ["var1", "var2"])
test_aux = AuxGradient(params)

# test aux
params = TestAuxParameters()
params.set("var", "testaux")
params.set("other_vars", ["var1", "var2"])
params.set("coefs", [1.5, 2.5])
aux = TestAux(params)

other_aux = [aux]
root_vars = ["grad_var1", "grad_var2"]
constant_data = {"var1": 1.2, "var2": 1.4}


class AuxGradientDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
