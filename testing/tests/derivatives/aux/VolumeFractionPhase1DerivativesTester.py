import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()
test_aux = factory.createObject("VolumeFractionPhase1", {"phase": 0})

other_aux = list()
root_vars = ["aA1"]


class VolumeFractionPhase1DerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
