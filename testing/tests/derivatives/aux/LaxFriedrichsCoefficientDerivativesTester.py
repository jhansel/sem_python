import unittest

from sem_python.aux.LaxFriedrichsCoefficient import LaxFriedrichsCoefficient, LaxFriedrichsCoefficientParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# LaxFriedrichs coefficient aux
params = LaxFriedrichsCoefficientParameters()
params.set("phase", 0)
params.set("var", "arhoA")
test_aux = LaxFriedrichsCoefficient(params)

# constant data
constant_data = {"dx": 0.5}

# positive velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arhoA1", "arhouA1"])
params.set("coefs", [2.0, 3.0])
u_positive_aux = TestAux(params)

# negative velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arhoA1", "arhouA1"])
params.set("coefs", [-1.2, -2.2])
u_negative_aux = TestAux(params)

# sound speed aux
params = TestAuxParameters()
params.set("var", "c1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
c_aux = TestAux(params)

other_aux_positive = [u_positive_aux, c_aux]
other_aux_negative = [u_negative_aux, c_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]

class LaxFriedrichsCoefficientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def testPositiveVelocity(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, other_aux_positive, root_vars, constant_data)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def testNegativeVelocity(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, other_aux_negative, root_vars, constant_data)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, other_aux_negative, root_vars, constant_data)
