import unittest

from sem_python.aux.SpecificVolume import SpecificVolume, SpecificVolumeParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# specific volume aux
params = SpecificVolumeParameters()
params.set("phase", 0)
test_aux = SpecificVolume(params)

# density aux
params = TestAuxParameters()
params.set("var", "rho1")
params.set("other_vars", ["aA1", "arhoA1"])
params.set("coefs", [2.0, 3.0])
rho_aux = TestAux(params)

other_aux = [rho_aux]
root_vars = ["aA1", "arhoA1"]

class SpecificVolumeDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
