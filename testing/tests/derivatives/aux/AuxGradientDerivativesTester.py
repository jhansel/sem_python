import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# gradient of test aux
params = dict()
params["aux"] = "testaux"
params["variable_names"] = ["var1", "var2"]
test_aux = factory.createObject("AuxGradient", params)

# test aux
params = dict()
params["var"] = "testaux"
params["other_vars"] = ["var1", "var2"]
params["coefs"] = [1.5, 2.5]
aux = factory.createObject("TestAux", params)

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
