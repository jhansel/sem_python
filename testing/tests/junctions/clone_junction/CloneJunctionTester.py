import unittest

import sem
from InputFileModification import InputFileModification
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class CloneJunctionTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "tests/junctions/clone_junction/"
    mods = list()
    mods.append(InputFileModification("NonlinearSolver", "verbose", False))
    mods.append(InputFileModification("Executioner", "verbose", False))
    mods.append(InputFileModification("Executioner", "end_time", 0.05))
    mods.append(InputFileModification("Output", "solution_file", "tests/junctions/clone_junction/solution_with_junction.csv"))
    mods.append(InputFileModification("Output", "plot_solution", False))
    sem.run(test_dir + "clone_junction_with_junction.in", mods)

    csv_tester = CSVTester(test_dir, "solution_with_junction.csv")
    self.assertTrue(csv_tester.filesAreEqual())

  def runDerivativeTest(self, test_weak):
    tester = JunctionTester("CloneJunction")
    rel_diffs = tester.checkJacobian(test_weak)
    n_i, n_j = rel_diffs.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)

  def testJacobianWeak(self):
    self.runDerivativeTest(True)

  def testJacobianStrong(self):
    self.runDerivativeTest(False)

if __name__ == "__main__":
  tester = JunctionTester("CloneJunction", verbose=True)
  _ = tester.checkJacobian(True)
  _ = tester.checkJacobian(False)
