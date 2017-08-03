import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from BerryInterfacialAreaDensity import BerryInterfacialAreaDensity, BerryInterfacialAreaDensityParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryInterfacialAreaDensityParameters()
params.set("a_int_min", 0.1)
params.set("a_int_max", 0.7)
test_aux = BerryInterfacialAreaDensity(params)
test_var = "a_int"

other_aux = dict()
other_vars = list()
root_vars = ["vf1"]

class BerryInterfacialAreaDensityDerivativesTester(unittest.TestCase):
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
