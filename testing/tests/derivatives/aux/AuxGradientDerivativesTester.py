import unittest

from AuxGradient import AuxGradient, AuxGradientParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# gradient of test aux
params = AuxGradientParameters()
params.set("aux", "testaux")
test_aux = AuxGradient(params)
test_var = "grad_testaux"

# test aux
params = TestAuxParameters()
params.set("var", "testaux")
params.set("other_vars", ["var1", "var2"])
params.set("coefs", [1.5, 2.5])
aux = TestAux(params)

other_aux = {"testaux" : aux}
other_vars = ["testaux"]
root_vars = ["grad_var1", "grad_var2"]
constant_data = {"var1": 1.2, "var2": 1.4}

class AuxGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux, other_vars, root_vars, constant_data)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux, other_vars, root_vars, constant_data)
