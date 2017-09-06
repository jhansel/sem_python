import unittest

import sem
from ParameterModification import BlockParameterModification, SubblockParameterModification
from CSVTester import CSVTester
from JunctionTester import JunctionTester

class EqualSolutionLM1PhaseJunctionTester(unittest.TestCase):
  def testSolution(self):
    input_dir = "tests/junctions/"
    test_dir = input_dir + "equal_solution_lm_1phase_junction/"
    solution_file = test_dir + "equal_solution_lm_1phase_junction.csv"
    mods = list()
    mods.append(SubblockParameterModification("Junctions", "junction1", "type", "EqualSolutionLM1PhaseJunction"))
    mods.append(SubblockParameterModification("Junctions", "junction1", "phase", "air"))
    mods.append(BlockParameterModification("NonlinearSolver", "verbose", False))
    mods.append(BlockParameterModification("Executioner", "verbose", False))
    mods.append(BlockParameterModification("Executioner", "end_time", 0.05))
    mods.append(BlockParameterModification("Output", "solution_file", solution_file))
    mods.append(BlockParameterModification("Output", "plot_solution", False))
    sem.run(input_dir + "junction.in", mods)

    csv_tester = CSVTester(test_dir, "equal_solution_lm_1phase_junction.csv")
    self.assertTrue(csv_tester.filesAreEqual())

  def runDerivativeTest(self, test_option):
    tester = JunctionTester("EqualSolutionLM1PhaseJunction")
    matched = tester.checkJacobian(test_option)
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])

  def testJacobianWeak(self):
    self.runDerivativeTest("weak")

  def testJacobianStrong(self):
    self.runDerivativeTest("strong")

if __name__ == "__main__":
  tester = JunctionTester("EqualSolutionLM1PhaseJunction", verbose=True)
  _ = tester.checkJacobian("weak")
  _ = tester.checkJacobian("strong")
