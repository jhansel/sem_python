import unittest

from sem_python.base.enums import ModelType
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

aux = {"g": list()}

class EnergyGravityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test1Phase(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyGravity", ModelType.OnePhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2PhaseNonInteracting(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyGravity", ModelType.TwoPhaseNonInteracting, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2Phase(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyGravity", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
