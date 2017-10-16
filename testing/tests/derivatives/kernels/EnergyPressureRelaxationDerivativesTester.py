import unittest

from sem_python.base.enums import ModelType
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

class EnergyPressureRelaxationDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()
    self.aux = {"pI_bar": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "p_relax": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"],
      "p1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"],
      "p2": ["aA1", "arhoA2", "arhouA2", "arhoEA2"]}

  def testPhase1(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 0, self.aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

  def testPhase2(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 1, self.aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)
