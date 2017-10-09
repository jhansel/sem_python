import unittest

from LaxFriedrichsCoefficientVolumeFraction import LaxFriedrichsCoefficientVolumeFraction, LaxFriedrichsCoefficientVolumeFractionParameters
from TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# Lax-Friedrichs coefficient aux
params = LaxFriedrichsCoefficientVolumeFractionParameters()
test_aux = LaxFriedrichsCoefficientVolumeFraction(params)

# constant data
constant_data = {"dx": 0.5}

# positive interfacial velocity aux
params = TestAuxParameters()
params.set("var", "uI")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5])
uI_positive_aux = TestAux(params)

# negative interfacial velocity aux
params = TestAuxParameters()
params.set("var", "uI")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [-1.3, -1.5, -1.7, -1.9, -2.1, -2.3, -2.5])
uI_negative_aux = TestAux(params)

other_aux_positive = [uI_positive_aux]
other_aux_negative = [uI_negative_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]

class LaxFriedrichsCoefficientVolumeFractionDerivativesTester(unittest.TestCase):
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
