import unittest

from sem_python.base.enums import ModelType
from KernelDerivativesTester import KernelDerivativesTester

class VolumeFractionPressureRelaxationDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test2Phase(self):
    aux = {"p_relax": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"],
      "p1": ["vf1", "arho1", "arhou1", "arhoE1"],
      "p2": ["vf1", "arho2", "arhou2", "arhoE2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("VolumeFractionPressureRelaxation", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"p_relax": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"],
    "p1": ["vf1", "arho1", "arhou1", "arhoE1"],
    "p2": ["vf1", "arho2", "arhou2", "arhoE2"]}
  _ = derivatives_tester.checkDerivatives("VolumeFractionPressureRelaxation", ModelType.TwoPhase, 0, aux)
