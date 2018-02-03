import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

# original aux
params = dict()
params["var"] = "aux_original"
params["other_vars"] = ["var1", "var2"]
params["coefs"] = [1.5, 2.5]
original_aux = factory.createObject("TestAux", params)

# copy aux
params = dict()
params["original_aux"] = "aux_original"
params["copy_aux"] = "aux_copy"
test_aux = factory.createObject("IdenticalAux", params)

other_aux = [original_aux]
root_vars = ["var1", "var2"]


class IdenticalAuxDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
