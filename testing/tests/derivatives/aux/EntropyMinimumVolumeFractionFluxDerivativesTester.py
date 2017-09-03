import unittest

from EntropyMinimumVolumeFractionFlux import EntropyMinimumVolumeFractionFlux, EntropyMinimumVolumeFractionFluxParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# viscous flux aux
params = EntropyMinimumVolumeFractionFluxParameters()
params.set("phase", 0)
test_aux = EntropyMinimumVolumeFractionFlux(params)

# viscous coefficient aux
params = TestAuxParameters()
params.set("var", "visccoef_vf")
params.set("other_vars", ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5, 2.2, 3.2, 4.2])
visccoef_vf_aux = TestAux(params)

# volume fraction gradient aux
params = TestAuxParameters()
params.set("var", "grad_vf")
params.set("other_vars", ["grad_vf1"])
params.set("coefs", [1.0])
grad_vf_aux = TestAux(params)

other_aux = [visccoef_vf_aux, grad_vf_aux]
root_vars = ["vf1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2", "grad_vf1"]

class EntropyMinimumVolumeFractionFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
