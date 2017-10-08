import unittest

from TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

params = TestAuxParameters()
params.set("var", "test")
params.set("other_vars", ["A", "B", "C"])
params.set("coefs", [2.0, 3.0, 4.0])
test_aux = TestAux(params)

other_aux = list()
root_vars = ["A", "B", "C"]

class TestAuxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
