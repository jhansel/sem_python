import unittest

from LaxFriedrichsCoefficient import LaxFriedrichsCoefficient, LaxFriedrichsCoefficientParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# LaxFriedrichs coefficient aux
params = LaxFriedrichsCoefficientParameters()
params.set("phase", 0)
params.set("var", "arho")
test_aux = LaxFriedrichsCoefficient(params)
test_var = "visccoef_arho1"

# constant data
constant_data = {"dx": 0.5}

# positive velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [2.0, 3.0])
u_positive_aux = TestAux(params)

# negative velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [-1.2, -2.2])
u_negative_aux = TestAux(params)

# sound speed aux
params = TestAuxParameters()
params.set("var", "c1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
c_aux = TestAux(params)

other_aux_positive = {"u1" : u_positive_aux, "c1" : c_aux}
other_aux_negative = {"u1" : u_negative_aux, "c1" : c_aux}
other_vars = ["u1", "c1"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1"]

class LaxFriedrichsCoefficientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def testPositiveVelocity(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux_positive, other_vars, root_vars, constant_data)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def testNegativeVelocity(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux_negative, other_vars, root_vars, constant_data)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux_negative, other_vars, root_vars, constant_data)
