import unittest

from enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class EnergyVolumeFractionGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def testPhase1(self):
    aux = {"uI": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "pI": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

  def testPhase2(self):
    aux = {"uI": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "pI": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 1, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"uI": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
    "pI": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]}
  _ = derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 0, aux)
