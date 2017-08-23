import unittest

import sem
from InputFileModification import InputFileModification
from CSVTester import CSVTester

class AmbrosoInterfaceClosuresTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/closures/ambroso_interface_closures/"
    sem.run(test_dir + "ambroso_interface_closures.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(InputFileModification("NonlinearSolver", "verbose", True))
  mods.append(InputFileModification("Executioner", "verbose", True))
  mods.append(InputFileModification("Output", "save_solution", False))
  sem.run("ambroso_interface_closures.in", mods)
