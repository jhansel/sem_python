import unittest

import sem
from InputFileModification import InputFileModification
from CSVTester import CSVTester

class MultiMeshTester(unittest.TestCase):
  ## Tests the solution
  def testSolution(self):
    test_dir = "tests/mesh/multi_mesh/"
    sem.run(test_dir + "multi_mesh.in")
    csv_tester = CSVTester(test_dir, "solution.csv")
    self.assertTrue(csv_tester.filesAreEqual())

if __name__ == "__main__":
  mods = list()
  mods.append(InputFileModification("NonlinearSolver", "verbose", True))
  mods.append(InputFileModification("Output", "solution_file", "solution.csv"))
  mods.append(InputFileModification("Output", "plot_solution", True))
  mods.append(InputFileModification("Output", "plot_file", "solution.pdf"))
  sem.run("multi_mesh.in", mods)
