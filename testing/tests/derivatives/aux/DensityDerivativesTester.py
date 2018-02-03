import unittest

from sem_python.base.Factory import Factory
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

factory = Factory()

test_aux = factory.createObject("Density", {"phase": 0})

vf_aux = factory.createObject("VolumeFractionPhase1", {"phase": 0})

other_aux = [vf_aux]
root_vars = ["aA1", "arhoA1"]


class DensityDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data={
                "A": 0.3
            })
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
