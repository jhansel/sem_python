import unittest

from sem_python.base.enums import ModelType
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

aux = {"uI": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
  "pI": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
  "vf1": ["aA1"]}

class EnergyVolumeFractionGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def testPhase1(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 0, aux, aux_gradients=["vf1"])
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

  def testPhase2(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 1, aux, aux_gradients=["vf1"])
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 0, aux, aux_gradients=["vf1"])
