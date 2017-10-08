import unittest

from ....src.testers.SolutionTester import SolutionTester

class SourceSplittingTester(unittest.TestCase):
  def test1Phase(self):
    test_dir = "testing/tests/solution/source_splitting/"
    input_file = test_dir + "source_splitting_1phase.in"
    solution_tester = SolutionTester(test_dir, input_file, solution_file_name="solution_1phase.csv")
    self.assertTrue(solution_tester.solutionsAreEqual())
