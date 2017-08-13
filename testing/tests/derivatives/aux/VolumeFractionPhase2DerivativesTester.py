import unittest

from VolumeFractionPhase2 import VolumeFractionPhase2, VolumeFractionPhase2Parameters
from AuxDerivativesTester import AuxDerivativesTester

params = VolumeFractionPhase2Parameters()
params.set("phase", 1)
test_aux = VolumeFractionPhase2(params)
test_var = "vf2"
other_aux = dict()
other_vars = list()
root_vars = ["vf1"]

class VolumeFractionPhase2DerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux, other_vars, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux, other_vars, root_vars)
