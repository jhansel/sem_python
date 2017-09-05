import unittest

from BerryInterfacialAreaDensity import BerryInterfacialAreaDensity, BerryInterfacialAreaDensityParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryInterfacialAreaDensityParameters()
params.set("a_int_min", 0.1)
params.set("a_int_max", 0.7)
test_aux = BerryInterfacialAreaDensity(params)

# phase-1 volume fraction aux
params = TestAuxParameters()
params.set("var", "vf1")
params.set("other_vars", ["aA1"])
params.set("coefs", [1.5])
vf1_aux = TestAux(params)

other_aux = [vf1_aux]
root_vars = ["aA1"]

class BerryInterfacialAreaDensityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
