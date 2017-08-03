import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from AmbrosoTheta import AmbrosoTheta, AmbrosoThetaParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# theta aux
params = AmbrosoThetaParameters()
params.set("pressure_relaxation_time", 2.0)
test_aux = AmbrosoTheta(params)
test_var = "theta"

# phase-1 pressure aux
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
p1_aux = TestAux(params)

# phase-2 pressure aux
params = TestAuxParameters()
params.set("var", "p2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
p2_aux = TestAux(params)

other_aux = {"p1" : p1_aux, "p2" : p2_aux}
other_vars = ["p1", "p2"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class AmbrosoThetaDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux, other_vars, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux, other_vars, root_vars)
