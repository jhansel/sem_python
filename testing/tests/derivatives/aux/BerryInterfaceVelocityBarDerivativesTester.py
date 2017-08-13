import unittest

from BerryInterfaceVelocityBar import BerryInterfaceVelocityBar, BerryInterfaceVelocityBarParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryInterfaceVelocityBarParameters()
test_aux = BerryInterfaceVelocityBar(params)
test_var = "uI_bar"

# phase-1 velocity
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [1.1, 1.2])
u1_aux = TestAux(params)

# phase-2 velocity
params = TestAuxParameters()
params.set("var", "u2")
params.set("other_vars", ["arho2", "arhou2"])
params.set("coefs", [1.7, 1.9])
u2_aux = TestAux(params)

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

other_aux = {"u1": u1_aux, "u2": u2_aux, "z1": z1_aux, "z2": z2_aux}
other_vars = ["u1", "u2", "z1", "z2"]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class BerryInterfaceVelocityBarDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(
      test_aux, test_var, other_aux, other_vars, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-5)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(
    test_aux, test_var, other_aux, other_vars, root_vars)
