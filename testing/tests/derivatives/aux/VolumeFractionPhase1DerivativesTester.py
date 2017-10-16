import unittest

from sem_python.aux.VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
test_aux = VolumeFractionPhase1(params)

other_aux = list()
root_vars = ["aA1"]

class VolumeFractionPhase1DerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
