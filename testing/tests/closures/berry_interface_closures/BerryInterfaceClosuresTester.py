import unittest

from ....src.testers.SolutionTester import SolutionTester


class BerryInterfaceClosuresTester(unittest.TestCase):

    def testSolution(self):
        test_dir = "testing/tests/closures/berry_interface_closures/"
        input_file = test_dir + "berry_interface_closures.in"
        solution_tester = SolutionTester(test_dir, input_file)
        self.assertTrue(solution_tester.solutionsAreEqual())
