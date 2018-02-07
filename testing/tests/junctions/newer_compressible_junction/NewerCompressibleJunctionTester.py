import unittest

from sem_python.input.InputFileModifier import InputFileModifier
from ....src.testers.SolutionTester import SolutionTester
from ....src.testers.JunctionTester import JunctionTester


class NewerCompressibleJunctionTester(unittest.TestCase):

    def testSolution(self):
        test_dir = "testing/tests/junctions/newer_compressible_junction/"
        input_file = "problems/junction_pressure_drop/junction_pressure_drop.in"

        input_file_modifier = InputFileModifier()
        input_file_modifier.modifySubblockParam("Executioner", "TimeStepSizer", "end_time", 0.05)

        solution_tester = SolutionTester(test_dir, input_file, input_file_modifier)
        self.assertTrue(solution_tester.solutionsAreEqual())

    def runDerivativeTest(self, test_option):
        tester = JunctionTester("NewerCompressibleJunction", rel_tol=1e-5)
        matched = tester.checkJacobian(
            test_option, junction_params={
                "loss_coefficients": [0.3, 0.4]
            })
        n_i, n_j = matched.shape
        for i in range(n_i):
            for j in range(n_j):
                self.assertTrue(matched[i, j])

    def testJacobianWeak(self):
        self.runDerivativeTest("weak")

    def testJacobianBoth(self):
        self.runDerivativeTest("both")
