import unittest

from VolumeFractionPhase2 import VolumeFractionPhase2, VolumeFractionPhase2Parameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = VolumeFractionPhase2Parameters()
params.set("phase", 1)
test_aux = VolumeFractionPhase2(params)

other_aux = list()
root_vars = ["aA1"]

class VolumeFractionPhase2DerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
