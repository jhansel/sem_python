import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

import unittest

sys.path.append(base_dir + "src/aux")
from BerryPressureRelaxationCoef import BerryPressureRelaxationCoef, BerryPressureRelaxationCoefParameters
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "testing/src/utilities")
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryPressureRelaxationCoefParameters()
test_aux = BerryPressureRelaxationCoef(params)
test_var = "p_relax"

# interfacial area density aux
params = TestAuxParameters()
params.set("var", "a_int")
params.set("other_vars", ["vf1"])
params.set("coefs", [1.5])
a_int_aux = TestAux(params)

# phase-1 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.2, 2.3, 3.4, 4.5])
z1_aux = TestAux(params)

# phase-2 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.1, 2.2, 3.3, 4.4])
z2_aux = TestAux(params)

other_aux = {"a_int": a_int_aux, "z1": z1_aux, "z2": z2_aux}
other_vars = ["a_int", "z1", "z2"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class BerryPressureRelaxationCoefDerivativesTester(unittest.TestCase):
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
