import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from AuxGradient import AuxGradient, AuxGradientParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# gradient of test aux
params = AuxGradientParameters()
params.set("aux", "testaux")
test_aux = AuxGradient(params)
test_var = "grad_testaux"

# test aux
params = TestAuxParameters()
params.set("var", "testaux")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [1.5, 2.5])
aux = TestAux(params)

other_aux = {"testaux" : aux}
other_vars = ["testaux"]
root_vars = ["arho1", "arhou1"]

class AuxGradientDerivativesTester(unittest.TestCase):
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
