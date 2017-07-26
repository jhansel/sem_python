import unittest

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "testing/src/utilities")
from KernelDerivativesTester import KernelDerivativesTester

class EnergyVolumeFractionGradientDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def testPhase1(self):
    aux = {"uI": ["arho1", "arhou1", "arho2", "arhou2"],
      "pI": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

  def testPhase2(self):
    aux = {"uI": ["arho1", "arhou1", "arho2", "arhou2"],
      "pI": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 1, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"uI": ["arho1", "arhou1", "arho2", "arhou2"],
    "pI": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"]}
  _ = derivatives_tester.checkDerivatives("EnergyVolumeFractionGradient", ModelType.TwoPhase, 0, aux)
