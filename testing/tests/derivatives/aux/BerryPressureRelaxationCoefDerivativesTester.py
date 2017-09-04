import unittest

from sem_python.aux.BerryPressureRelaxationCoef import BerryPressureRelaxationCoef, BerryPressureRelaxationCoefParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = BerryPressureRelaxationCoefParameters()
test_aux = BerryPressureRelaxationCoef(params)

# interfacial area density aux
params = TestAuxParameters()
params.set("var", "a_int")
params.set("other_vars", ["vf1"])
params.set("coefs", [1.5])
a_int_aux = TestAux(params)

# phase-1 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z1")
params.set("other_vars", ["vf1", "arho1", "arhou1", "arhoE1"])
params.set("coefs", [1.2, 2.3, 3.4, 4.5])
z1_aux = TestAux(params)

# phase-2 acoustic impedance aux
params = TestAuxParameters()
params.set("var", "z2")
params.set("other_vars", ["vf1", "arho2", "arhou2", "arhoE2"])
params.set("coefs", [1.1, 2.2, 3.3, 4.4])
z2_aux = TestAux(params)

other_aux = [a_int_aux, z1_aux, z2_aux]
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class BerryPressureRelaxationCoefDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
