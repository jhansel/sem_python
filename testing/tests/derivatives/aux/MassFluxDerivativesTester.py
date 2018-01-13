import unittest

from sem_python.aux.MassFlux import MassFlux, MassFluxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = MassFluxParameters()
params.set("phase", 0)
test_aux = MassFlux(params)

other_aux = list()
root_vars = ["arhouA1"]


class MassFluxDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            test_aux, other_aux, root_vars, constant_data={
                "A": 0.3
            })
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
