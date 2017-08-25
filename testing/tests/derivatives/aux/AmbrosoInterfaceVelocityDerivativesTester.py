import unittest

from AmbrosoInterfaceVelocity import AmbrosoInterfaceVelocity, AmbrosoInterfaceVelocityParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# interface velocity aux
params = AmbrosoInterfaceVelocityParameters()
test_aux = AmbrosoInterfaceVelocity(params)

# phase-1 velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [1.2, 2.2])
u1_aux = TestAux(params)

# phase-2 velocity aux
params = TestAuxParameters()
params.set("var", "u2")
params.set("other_vars", ["arho2", "arhou2"])
params.set("coefs", [1.5, 2.5])
u2_aux = TestAux(params)

# beta aux
params = TestAuxParameters()
params.set("var", "beta")
params.set("other_vars", ["arho1", "arho2"])
params.set("coefs", [1.7, 3.2])
beta_aux = TestAux(params)

other_aux = [u1_aux, u2_aux, beta_aux]
root_vars = ["arho1", "arhou1", "arho2", "arhou2"]

class AmbrosoInterfaceVelocityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
