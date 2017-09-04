import unittest

from sem_python.aux.ConstantAux import ConstantAux, ConstantAuxParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# test aux
params = ConstantAuxParameters()
params.set("name", "test_aux")
params.set("value", 6.0)
test_aux = ConstantAux(params)

other_aux = list()
root_vars = ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]

class ConstantAuxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
