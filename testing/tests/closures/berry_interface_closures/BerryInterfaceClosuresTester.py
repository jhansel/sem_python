import unittest

import sem
from ParameterModification import BlockParameterModification
from CSVTester import CSVTester

class BerryInterfaceClosuresTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/closures/berry_interface_closures/"
    sem.run(test_dir + "berry_interface_closures.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(BlockParameterModification("NonlinearSolver", "verbose", True))
  mods.append(BlockParameterModification("Executioner", "verbose", True))
  mods.append(BlockParameterModification("Output", "save_solution", False))
  sem.run("berry_interface_closures.in", mods)
