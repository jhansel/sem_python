import unittest

from ...src.testers.BCTester import BCTester

class FreeBCTester(unittest.TestCase):
  def testJacobian(self):
    self.tester = BCTester()
    rel_diffs = self.tester.checkJacobian("FreeBC", bc_params={"phase": "0"})
    n_i, n_j = rel_diffs.shape
    for i in range(n_i):
      for j in range(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)
