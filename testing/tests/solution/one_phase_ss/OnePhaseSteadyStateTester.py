import unittest

from SolutionTester import SolutionTester

class OnePhaseSteadyStateTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "tests/solution/one_phase_ss/"
    input_file = test_dir + "one_phase_ss.in"
    solution_tester = SolutionTester(test_dir, input_file)
    self.assertTrue(solution_tester.solutionsAreEqual())
