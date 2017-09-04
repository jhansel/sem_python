import unittest

from sem_python.aux.BerryInterfaceVelocity import BerryInterfaceVelocity, BerryInterfaceVelocityParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryInterfaceVelocityParameters()
test_aux = BerryInterfaceVelocity(params)

# bar interface velocity aux
params = TestAuxParameters()
params.set("var", "uI_bar")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
uI_bar_aux = TestAux(params)

# phase-1 pressure aux
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [2.1, 2.3, 3.4, 4.5])
p1_aux = TestAux(params)

# phase-2 pressure aux
params = TestAuxParameters()
params.set("var", "p2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [2.4, 2.2, 3.3, 4.4])
p2_aux = TestAux(params)

# phase-1 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.6, 2.3, 4.5, 2.1])
z1_aux = TestAux(params)

# phase-2 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.2, 3.2, 4.1, 2.4])
z2_aux = TestAux(params)

other_aux = [uI_bar_aux, p1_aux, p2_aux, z1_aux, z2_aux]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]
constant_data_positive = {"grad_vf1": 0.6}
constant_data_negative = {"grad_vf1": -0.6}

class BerryInterfaceVelocityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def testPositiveGradient(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data_positive)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def testNegativeGradient(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data_negative)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data_negative)
