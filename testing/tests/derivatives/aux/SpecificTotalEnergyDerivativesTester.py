import unittest

from sem_python.aux.SpecificTotalEnergy import SpecificTotalEnergy, SpecificTotalEnergyParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = SpecificTotalEnergyParameters()
params.set("phase", 0)
test_aux = SpecificTotalEnergy(params)

other_aux = list()
root_vars = ["arhoA1", "arhoEA1"]


class SpecificTotalEnergyDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = AuxDerivativesTester()

    def test(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
