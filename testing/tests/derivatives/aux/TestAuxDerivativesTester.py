import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()
params = dict()
params["var"] = "test"
params["other_vars"] = ["A", "B", "C"]
params["coefs"] = [2.0, 3.0, 4.0]
test_aux = factory.createObject("TestAux", params)

other_aux = list()
root_vars = ["A", "B", "C"]


class TestAuxDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
