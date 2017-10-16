import unittest

from sem_python.aux.BerryInterfacePressureBar import BerryInterfacePressureBar, BerryInterfacePressureBarParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryInterfacePressureBarParameters()
test_aux = BerryInterfacePressureBar(params)

# phase-1 pressure
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.1, 1.2, 1.3, 1.4])
p1_aux = TestAux(params)

# phase-2 pressure
params = TestAuxParameters()
params.set("var", "p2")
params.set("other_vars", ["aA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.7, 1.9, 2.3, 1.6])
p2_aux = TestAux(params)

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

other_aux = [p1_aux, p2_aux, z1_aux, z2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]

class BerryInterfacePressureBarDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)
