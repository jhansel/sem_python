import unittest

from sem_python.aux.InternalEnergyDensity import InternalEnergyDensity, InternalEnergyDensityParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# internal energy density aux
params = InternalEnergyDensityParameters()
params.set("phase", 0)
test_aux = InternalEnergyDensity(params)

# density aux
params = TestAuxParameters()
params.set("var", "rho1")
params.set("other_vars", ["aA1", "arhoA1"])
params.set("coefs", [2.0, 3.0])
rho_aux = TestAux(params)

# specific internal energy aux
params = TestAuxParameters()
params.set("var", "e1")
params.set("other_vars", ["arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [2.5, 3.5, 4.5])
e_aux = TestAux(params)

other_aux = [rho_aux, e_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]

class InternalEnergyDensityDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
