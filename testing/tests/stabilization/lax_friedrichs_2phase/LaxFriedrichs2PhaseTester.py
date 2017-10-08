import unittest

from ....src.testers.SolutionTester import SolutionTester

class LaxFriedrichs2PhaseTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "testing/tests/stabilization/lax_friedrichs_2phase/"
    input_file = test_dir + "lax_friedrichs_2phase.in"
    solution_tester = SolutionTester(test_dir, input_file)
    self.assertTrue(solution_tester.solutionsAreEqual())
