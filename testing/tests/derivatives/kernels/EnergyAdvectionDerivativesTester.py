import unittest

from enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class EnergyAdvectionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test1Phase(self):
    aux = {"u1": ["arho1", "arhou1"], "vf1": list(), "p1": ["arho1", "arhou1", "arhoE1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyAdvection", ModelType.OnePhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2PhaseNonInteracting(self):
    aux = {"u1": ["arho1", "arhou1"], "vf1": list(), "p1": ["arho1", "arhou1", "arhoE1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyAdvection", ModelType.TwoPhaseNonInteracting, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2Phase(self):
    aux = {"u1": ["arho1", "arhou1"], "vf1": ["vf1"], "p1": ["vf1", "arho1", "arhou1", "arhoE1"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyAdvection", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"u1": ["arho1", "arhou1"], "vf1": list(), "p1": ["arho1", "arhou1", "arhoE1"]}
  _ = derivatives_tester.checkDerivatives("EnergyAdvection", ModelType.TwoPhaseNonInteracting, 0, aux)
