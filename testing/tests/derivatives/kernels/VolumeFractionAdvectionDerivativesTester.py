import unittest

from enums import ModelType
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

class VolumeFractionAdvectionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test2Phase(self):
    aux = {"uI": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"], "vf1": ["aA1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("VolumeFractionAdvection", ModelType.TwoPhase, 0, aux, aux_gradients=["vf1"])
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"uI": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"], "vf1": ["aA1"]}
  _ = derivatives_tester.checkDerivatives("VolumeFractionAdvection", ModelType.TwoPhase, 0, aux, aux_gradients=["vf1"])
