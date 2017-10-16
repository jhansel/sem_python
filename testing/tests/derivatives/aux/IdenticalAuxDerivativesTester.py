import unittest

from sem_python.aux.IdenticalAux import IdenticalAux, IdenticalAuxParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# original aux
params = TestAuxParameters()
params.set("var", "aux_original")
params.set("other_vars", ["var1", "var2"])
params.set("coefs", [1.5, 2.5])
original_aux = TestAux(params)

# copy aux
params = IdenticalAuxParameters()
params.set("original_aux", "aux_original")
params.set("copy_aux", "aux_copy")
test_aux = IdenticalAux(params)

other_aux = [original_aux]
root_vars = ["var1", "var2"]

class IdenticalAuxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
