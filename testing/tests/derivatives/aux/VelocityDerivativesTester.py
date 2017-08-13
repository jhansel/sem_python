import unittest

from Velocity import Velocity, VelocityParameters
from AuxDerivativesTester import AuxDerivativesTester

params = VelocityParameters()
params.set("phase", 0)
test_aux = Velocity(params)
test_var = "u1"
other_aux = dict()
other_vars = list()
root_vars = ["arho1", "arhou1"]

class VelocityDerivativesTester(unittest.TestCase):
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
