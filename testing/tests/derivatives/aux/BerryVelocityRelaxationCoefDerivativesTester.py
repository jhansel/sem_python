import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from BerryVelocityRelaxationCoef import BerryVelocityRelaxationCoef, BerryVelocityRelaxationCoefParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryVelocityRelaxationCoefParameters()
test_aux = BerryVelocityRelaxationCoef(params)
test_var = "u_relax"

# pressure relaxation aux
params = TestAuxParameters()
params.set("var", "p_relax")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
p_relax_aux = TestAux(params)

# phase-1 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [2.1, 2.3, 3.4, 4.5])
z1_aux = TestAux(params)

# phase-2 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [2.4, 2.2, 3.3, 4.4])
z2_aux = TestAux(params)

other_aux = {"p_relax": p_relax_aux, "z1": z1_aux, "z2": z2_aux}
other_vars = ["p_relax", "z1", "z2"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class BerryVelocityRelaxationCoefDerivativesTester(unittest.TestCase):
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
