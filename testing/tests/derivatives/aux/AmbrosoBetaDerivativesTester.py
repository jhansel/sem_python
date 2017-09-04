import unittest

from sem_python.aux.AmbrosoBeta import AmbrosoBeta, AmbrosoBetaParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# beta aux
params = AmbrosoBetaParameters()
params.set("chi", 0.5)
test_aux = AmbrosoBeta(params)

other_aux = []
root_vars = ["arho1", "arho2"]

class AmbrosoBetaDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
