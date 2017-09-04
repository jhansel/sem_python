import unittest

from sem_python.aux.SpecificInternalEnergy import SpecificInternalEnergy, SpecificInternalEnergyParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# specific internal energy aux
params = SpecificInternalEnergyParameters()
params.set("phase", 0)
test_aux = SpecificInternalEnergy(params)

# velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arho1", "arhou1"])
params.set("coefs", [2.0, 3.0])
u_aux = TestAux(params)

# specific total energy aux
params = TestAuxParameters()
params.set("var", "E1")
params.set("other_vars", ["arho1", "arhoE1"])
params.set("coefs", [2.5, 3.5])
E_aux = TestAux(params)

other_aux = [u_aux, E_aux]
root_vars = ["arho1", "arhou1", "arhoE1"]

class SpecificInternalEnergyDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
