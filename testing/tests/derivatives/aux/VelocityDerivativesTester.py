import unittest

from sem_python.aux.Velocity import Velocity, VelocityParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = VelocityParameters()
params.set("phase", 0)
test_aux = Velocity(params)

other_aux = list()
root_vars = ["arhoA1", "arhouA1"]

class VelocityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
