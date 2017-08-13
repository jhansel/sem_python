import unittest

from Density import Density, DensityParameters
from VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from AuxDerivativesTester import AuxDerivativesTester

params = DensityParameters()
params.set("phase", 0)
test_aux = Density(params)

test_var = "rho1"

params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)
other_aux = {"vf1" : vf_aux}

other_vars = ["vf1"]

root_vars = ["vf1", "arho1"]

class DensityDerivativesTester(unittest.TestCase):
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
