import unittest

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "testing/src/utilities")
from KernelDerivativesTester import KernelDerivativesTester

class VolumeFractionAdvectionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()

  def test2Phase(self):
    aux = {"uI": ["arho1", "arhou1", "arho2", "arhou2"]}
    rel_diffs = self.derivatives_tester.checkDerivatives("VolumeFractionAdvection", ModelType.TwoPhase, 0, aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"uI": ["arho1", "arhou1", "arho2", "arhou2"]}
  _ = derivatives_tester.checkDerivatives("VolumeFractionAdvection", ModelType.TwoPhase, 0, aux)
