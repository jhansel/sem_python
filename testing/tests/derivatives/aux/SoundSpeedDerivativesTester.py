import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from SoundSpeed import SoundSpeed, SoundSpeedParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

def computeSoundSpeed(v, e):
  v_slope = 2.0
  e_slope = 3.0
  c = v_slope * v + e_slope * e
  return (c, v_slope, e_slope)

# sound speed aux
params = SoundSpeedParameters()
params.set("phase", 0)
params.set("c_function", computeSoundSpeed)
test_aux = SoundSpeed(params)
test_var = "c1"

# specific volume aux
params = TestAuxParameters()
params.set("var", "v1")
params.set("other_vars", ["vf1", "arho1"])
params.set("coefs", [2.0, 3.0])
v_aux = TestAux(params)

# pressure aux
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
p_aux = TestAux(params)

other_aux = {"v1" : v_aux, "p1" : p_aux}
other_vars = ["v1", "p1"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1"]

class SoundSpeedDerivativesTester(unittest.TestCase):
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
