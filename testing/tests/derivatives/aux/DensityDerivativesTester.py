import unittest

from sem_python.aux.Density import Density, DensityParameters
from sem_python.aux.VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = DensityParameters()
params.set("phase", 0)
test_aux = Density(params)

params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)

other_aux = [vf_aux]
root_vars = ["aA1", "arhoA1"]

class DensityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data={"A": 0.3})
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
