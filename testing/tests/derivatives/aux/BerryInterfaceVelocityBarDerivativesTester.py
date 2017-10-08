import unittest

from BerryInterfaceVelocityBar import BerryInterfaceVelocityBar, BerryInterfaceVelocityBarParameters
from TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryInterfaceVelocityBarParameters()
test_aux = BerryInterfaceVelocityBar(params)

# phase-1 velocity
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arhoA1", "arhouA1"])
params.set("coefs", [1.1, 1.2])
u1_aux = TestAux(params)

# phase-2 velocity
params = TestAuxParameters()
params.set("var", "u2")
params.set("other_vars", ["arhoA2", "arhouA2"])
params.set("coefs", [1.7, 1.9])
u2_aux = TestAux(params)

# phase-1 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [2.1, 2.3, 3.4, 4.5])
z1_aux = TestAux(params)

# phase-2 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z2")
params.set("other_vars", ["aA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [2.4, 2.2, 3.3, 4.4])
z2_aux = TestAux(params)

other_aux = [u1_aux, u2_aux, z1_aux, z2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]

class BerryInterfaceVelocityBarDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-5)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
