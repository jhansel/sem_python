import unittest

import sem
from ParameterModification import BlockParameterModification
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
  mods.append(BlockParameterModification("NonlinearSolver", "verbose", True))
  mods.append(BlockParameterModification("Output", "solution_file", "solution.csv"))
  mods.append(BlockParameterModification("Output", "plot_solution", True))
  mods.append(BlockParameterModification("Output", "plot_file", "solution.pdf"))
  sem.run("multi_mesh.in", mods)
