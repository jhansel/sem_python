import unittest

from sem_python.aux.EntropyMinimumVolumeFractionFlux import EntropyMinimumVolumeFractionFlux, EntropyMinimumVolumeFractionFluxParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumVolumeFractionFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumVolumeFractionFlux(params)

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_aA1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5, 2.2, 3.2, 4.2])
visccoef_aA1_aux = TestAux(params)

# volume fraction gradient aux
params = TestAuxParameters()
params.set("var", "grad_vf1")
params.set("other_vars", ["grad_aA1"])
params.set("coefs", [1.0])
grad_vf_aux = TestAux(params)

other_aux = [visccoef_aA1_aux, grad_vf_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_aA1"]

class EntropyMinimumVolumeFractionFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data={"A": 0.3})
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data={"A": 0.3})
