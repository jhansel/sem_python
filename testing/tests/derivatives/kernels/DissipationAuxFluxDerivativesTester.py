import unittest

from enums import ModelType, VariableName
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

class DissipationAuxFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test1Phase(self):
    aux = {"viscflux_arhouA1": ["arhoA1", "arhouA1", "arhoEA1", "grad_arhoA1", "grad_arhouA1", "grad_arhoEA1"]}
    kernel_params = {"var_enum": VariableName.ARhoUA, "flux_name": "viscflux_arhouA1"}
    rel_diffs = self.derivatives_tester.checkDerivatives("DissipationAuxFlux", ModelType.OnePhase, 0, aux, kernel_params=kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2Phase(self):
    aux = {"viscflux_arhouA1": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "grad_aA1", "grad_arhoA1", "grad_arhouA1", "grad_arhoEA1"]}
    kernel_params = {"var_enum": VariableName.ARhoUA, "flux_name": "viscflux_arhouA1"}
    rel_diffs = self.derivatives_tester.checkDerivatives("DissipationAuxFlux", ModelType.TwoPhase, 0, aux, kernel_params=kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  aux = {"viscflux_arhouA1": ["aA1", "arhoA1", "arhouA1", "arhoEA1", "grad_aA1", "grad_arhoA1", "grad_arhouA1", "grad_arhoEA1"]}
  kernel_params = {"var_enum": VariableName.ARhoUA, "flux_name": "viscflux_arhouA1"}
  derivatives_tester = KernelDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives("DissipationAuxFlux", ModelType.TwoPhase, 0, aux, kernel_params=kernel_params)
