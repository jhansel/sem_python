import unittest

from sem_python.base.enums import ModelType, VariableName
from ....src.testers.KernelDerivativesTester import KernelDerivativesTester

class ArbitraryAuxKernel1PhaseDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test(self):
    aux = {"var1": ["arhoA1", "arhouA1"]}
    kernel_params = {"var_enum": VariableName.ARhoA, "aux_name": "var1"}
    rel_diffs = self.derivatives_tester.checkDerivatives("ArbitraryAuxKernel1Phase", \
      ModelType.OnePhase, 0, aux, kernel_params=kernel_params)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
