import unittest

from sem_python.base.enums import ModelType
from ...src.testers.BCTester import BCTester


class DirichletVolumeFractionBCTester(unittest.TestCase):

    def setUp(self):
        self.tester = BCTester()

    def testJacobian(self):
        rel_diffs = self.tester.checkJacobian(
            "DirichletVolumeFractionBC", model_type=ModelType.TwoPhase, bc_params={
                "vf1": 0.1
            })
        n_i, n_j = rel_diffs.shape
        for i in range(n_i):
            for j in range(n_j):
                self.assertLessEqual(rel_diffs[i, j], 1e-6)
