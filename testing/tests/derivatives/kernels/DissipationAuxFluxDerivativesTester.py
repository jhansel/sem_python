import unittest

from sem_python.base.enums import ModelType, VariableName
from KernelDerivativesTester import KernelDerivativesTester

class DissipationAuxFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test1Phase(self):
    aux = {"viscflux_arhou1": ["arho1", "arhou1", "arhoE1", "grad_arho1", "grad_arhou1", "grad_arhoE1"]}
    kernel_params = {"var_enum": VariableName.ARhoU, "flux_name": "viscflux_arhou1"}
    rel_diffs = self.derivatives_tester.checkDerivatives("DissipationAuxFlux", ModelType.OnePhase, 0, aux, kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2Phase(self):
    aux = {"viscflux_arhou1": ["vf1", "arho1", "arhou1", "arhoE1", "grad_vf1", "grad_arho1", "grad_arhou1", "grad_arhoE1"]}
    kernel_params = {"var_enum": VariableName.ARhoU, "flux_name": "viscflux_arhou1"}
    rel_diffs = self.derivatives_tester.checkDerivatives("DissipationAuxFlux", ModelType.TwoPhase, 0, aux, kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

if __name__ == "__main__":
  aux = {"viscflux_arhou1": ["vf1", "arho1", "arhou1", "arhoE1", "grad_vf1", "grad_arho1", "grad_arhou1", "grad_arhoE1"]}
  kernel_params = {"var_enum": VariableName.ARhoU, "flux_name": "viscflux_arhou1"}
  derivatives_tester = KernelDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives("DissipationAuxFlux", ModelType.TwoPhase, 0, aux, kernel_params)
