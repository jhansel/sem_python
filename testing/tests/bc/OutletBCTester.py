import unittest

from ...src.testers.BCTester import BCTester

class OutletBCTester(unittest.TestCase):
  def setUp(self):
    self.tester = BCTester()

  def testJacobianWeaklyEnforce(self):
    rel_diffs = self.tester.checkJacobian("OutletBC", bc_params={"phase": "0", "p": "1.2",
      "strongly_enforce_energy": False})
    n_i, n_j = rel_diffs.shape
    for i in range(n_i):
      for j in range(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)

  def testJacobianStronglyEnforce(self):
    rel_diffs = self.tester.checkJacobian("OutletBC", bc_params={"phase": "0", "p": "1.2",
      "strongly_enforce_energy": True})
    n_i, n_j = rel_diffs.shape
    for i in range(n_i):
      for j in range(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)
