import unittest

import sem
from sem_python.input.InputFileModifier import InputFileModifier
from ....src.testers.CSVTester import CSVTester
from ....src.testers.JunctionTester import JunctionTester


class NewCompressibleJunctionTester(unittest.TestCase):

    def runDerivativeTest(self, test_option):
        tester = JunctionTester("NewCompressibleJunction")
        matched = tester.checkJacobian(test_option)
        n_i, n_j = matched.shape
        for i in range(n_i):
            for j in range(n_j):
                self.assertTrue(matched[i, j])

    def testJacobianWeak(self):
        self.runDerivativeTest("weak")

    def testJacobianBoth(self):
        self.runDerivativeTest("both")
