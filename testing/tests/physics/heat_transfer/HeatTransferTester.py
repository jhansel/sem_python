import unittest

from ....src.testers.SolutionTester import SolutionTester

class HeatTransferTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "testing/tests/physics/heat_transfer/"
    input_file = "problems/heat_transfer/heat_transfer.in"

    solution_tester = SolutionTester(test_dir, input_file)
    self.assertTrue(solution_tester.solutionsAreEqual())
