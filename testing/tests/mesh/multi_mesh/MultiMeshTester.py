import unittest

from SolutionTester import SolutionTester

class MultiMeshTester(unittest.TestCase):
  def testSolution(self):
    test_dir = "tests/mesh/multi_mesh/"
    input_file = test_dir + "multi_mesh.in"
    solution_tester = SolutionTester(test_dir, input_file)
    self.assertTrue(solution_tester.solutionsAreEqual())
