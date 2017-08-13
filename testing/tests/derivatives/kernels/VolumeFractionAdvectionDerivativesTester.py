import unittest

from enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class VolumeFractionAdvectionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test2Phase(self):
    aux = {"uI": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("VolumeFractionAdvection", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"uI": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]}
  _ = derivatives_tester.checkDerivatives("VolumeFractionAdvection", ModelType.TwoPhase, 0, aux)
