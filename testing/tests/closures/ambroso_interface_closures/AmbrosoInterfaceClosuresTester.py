import unittest
import filecmp

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir)
import sem

sys.path.append(base_dir + "src/input")
from InputFileModification import InputFileModification

class AmbrosoInterfaceClosuresTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/closures/ambroso_interface_closures/"
    sem.run(test_dir + "ambroso_interface_closures.in")
    self.assertTrue(filecmp.cmp(test_dir + "solution.csv", test_dir + "gold/solution.csv"))

if __name__ == "__main__":
  mods = list()
  mods.append(InputFileModification("NonlinearSolver", "verbose", True))
  mods.append(InputFileModification("Executioner", "verbose", True))
  mods.append(InputFileModification("Output", "save_solution", False))
  sem.run("ambroso_interface_closures.in", mods)
