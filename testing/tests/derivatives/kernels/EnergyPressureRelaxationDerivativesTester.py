import unittest

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "testing/src/utilities")
from KernelDerivativesTester import KernelDerivativesTester

class EnergyPressureRelaxationDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = KernelDerivativesTester()
    self.aux = {"pI_bar": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"],
      "p_relax": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"],
      "p1": ["vf1", "arho1", "arhou1", "arhoE1"],
      "p2": ["vf1", "arho2", "arhou2", "arhoE2"]}

  def testPhase1(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 0, self.aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

  def testPhase2(self):
    rel_diffs = self.derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 1, self.aux)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-5)

if __name__ == "__main__":
  derivatives_tester = KernelDerivativesTester(True)
  aux = {"pI_bar": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"],
    "p_relax": ["vf1", "arho1", "arhou1", "arhoE1", "arho2", "arhou2", "arhoE2"],
    "p1": ["vf1", "arho1", "arhou1", "arhoE1"],
    "p2": ["vf1", "arho2", "arhou2", "arhoE2"]}
  _ = derivatives_tester.checkDerivatives("EnergyPressureRelaxation", ModelType.TwoPhase, 0, aux)
