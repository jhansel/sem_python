import unittest

from ...src.testers.BCTester import BCTester

bc_params={"phase": "0"}

class SolidWallBCTester(unittest.TestCase):
  def setUp(self):
    self.tester = BCTester()

  def testJacobian(self):
    rel_diffs = self.tester.checkJacobian("SolidWallBC", bc_params=bc_params)
    n_i, n_j = rel_diffs.shape
    for i in range(n_i):
      for j in range(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)
