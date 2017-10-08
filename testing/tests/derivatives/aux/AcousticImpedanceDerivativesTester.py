import unittest

from sem_python.aux.AcousticImpedance import AcousticImpedance, AcousticImpedanceParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# test aux
params = AcousticImpedanceParameters()
params.set("phase", 0)
test_aux = AcousticImpedance(params)

# density aux
params = TestAuxParameters()
params.set("var", "rho1")
params.set("other_vars", ["aA1", "arhoA1"])
params.set("coefs", [1.4, 2.5])
rho_aux = TestAux(params)

# sound speed aux
params = TestAuxParameters()
params.set("var", "c1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
c_aux = TestAux(params)

other_aux = [rho_aux, c_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]

class AcousticImpedanceDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
