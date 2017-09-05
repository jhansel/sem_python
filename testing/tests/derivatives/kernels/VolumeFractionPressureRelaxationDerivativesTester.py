import unittest

from enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class VolumeFractionPressureRelaxationDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test2Phase(self):
    aux = {"p_relax": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "p1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"],
      "p2": ["aA1", "arhoA2", "arhouA2", "arhoEA2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("VolumeFractionPressureRelaxation", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"p_relax": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
    "p1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"],
    "p2": ["aA1", "arhoA2", "arhouA2", "arhoEA2"]}
  _ = derivatives_tester.checkDerivatives("VolumeFractionPressureRelaxation", ModelType.TwoPhase, 0, aux)
