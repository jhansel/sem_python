import unittest

from enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class EnergyPressureRelaxationDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()
    self.aux = {"pI_bar": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "p_relax": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "p1": ["vf1", "arhoA1", "arhouA1", "arhoEA1"],
      "p2": ["vf1", "arhoA2", "arhouA2", "arhoEA2"]}

  def testPhase1(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 0, self.aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

  def testPhase2(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 1, self.aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"pI_bar": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
    "p_relax": ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
    "p1": ["vf1", "arhoA1", "arhouA1", "arhoEA1"],
    "p2": ["vf1", "arhoA2", "arhouA2", "arhoEA2"]}
  _ = derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 0, aux)
