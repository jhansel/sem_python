import unittest

from SolutionTester import SolutionTester

class TwoPhaseNoninteractingSteadyStateTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "tests/solution/two_phase_noninteracting_ss/"
    input_file = test_dir + "two_phase_noninteracting_ss.in"
    solution_tester = SolutionTester(test_dir, input_file)
    self.assertTrue(solution_tester.solutionsAreEqual())
