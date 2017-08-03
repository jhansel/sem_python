import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from AmbrosoBeta import AmbrosoBeta, AmbrosoBetaParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# beta aux
params = AmbrosoBetaParameters()
params.set("chi", 0.5)
test_aux = AmbrosoBeta(params)
test_var = "beta"

other_aux = {}
other_vars = []
root_vars = ["arho1", "arho2"]

class AmbrosoBetaDerivativesTester(unittest.TestCase):
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