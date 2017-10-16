import unittest

from sem_python.base.enums import ModelType
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

class MomentumVolumeFractionGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def testPhase1(self):
    aux = {"pI": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"], "vf1": ["aA1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("MomentumVolumeFractionGradient", ModelType.TwoPhase, 0, aux, aux_gradients=["vf1"])
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

  def testPhase2(self):
    aux = {"pI": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"], "vf1": ["aA1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("MomentumVolumeFractionGradient", ModelType.TwoPhase, 1, aux, aux_gradients=["vf1"])
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)
