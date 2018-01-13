import unittest

from sem_python.base.enums import ModelType
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

aux = {"T1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"]}


class EnergyHeatTransferDerivativesTester(unittest.TestCase):

    def setUp(self):
        self.derivatives_tester = KernelDerivativesTester()

    def test1Phase(self):
        rel_diffs = self.derivatives_tester.checkDerivatives(
            "EnergyHeatTransfer", ModelType.OnePhase, 0, aux)
        for key in rel_diffs:
            self.assertLessEqual(rel_diffs[key], 1e-6)
