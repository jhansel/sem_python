import unittest

from sem_python.base.enums import ModelType, VariableName
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

class DissipationVariableGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test1Phase(self):
    aux = {"visccoef_arhouA1": ["arhoA1", "arhouA1", "arhoEA1"]}
    kernel_params = {"var_enum": VariableName.ARhoUA}
    rel_diffs = self.derivatives_tester.checkDerivatives("DissipationVariableGradient", ModelType.OnePhase, 0, aux, kernel_params=kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)

  def test2Phase(self):
    aux = {"visccoef_aA1": ["aA1", "arhoA1", "arhouA1", "arhoEA1"]}
    kernel_params = {"var_enum": VariableName.AA1}
    rel_diffs = self.derivatives_tester.checkDerivatives("DissipationVariableGradient", ModelType.TwoPhase, 0, aux, kernel_params=kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)
