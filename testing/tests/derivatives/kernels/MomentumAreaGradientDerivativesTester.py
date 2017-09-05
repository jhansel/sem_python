import unittest

from enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class MomentumAreaGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test1Phase(self):
    aux = {"vf1": list(), "p1": ["arhoA1", "arhouA1", "arhoEA1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("MomentumAreaGradient", ModelType.OnePhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2PhaseNonInteracting(self):
    aux = {"vf1": list(), "p1": ["arhoA1", "arhouA1", "arhoEA1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("MomentumAreaGradient", ModelType.TwoPhaseNonInteracting, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2Phase(self):
    aux = {"vf1": ["aA1"], "p1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("MomentumAreaGradient", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"vf1": ["aA1"], "p1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"]}
  _ = derivatives_tester.checkDerivatives("MomentumAreaGradient", ModelType.TwoPhase, 0, aux)
