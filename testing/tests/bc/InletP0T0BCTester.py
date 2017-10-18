import unittest

from ...src.testers.BCTester import BCTester

class InletP0T0BCTester(unittest.TestCase):
  def setUp(self):
    self.tester = BCTester()

  def testJacobianWeaklyEnforce(self):
    rel_diffs = self.tester.checkJacobian("InletP0T0BC", bc_params={"phase": "0", "p0": "0.2", "T0": "0.6", "strongly_enforce_energy": False})
    n_i, n_j = rel_diffs.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)

  def testJacobianStronglyEnforce(self):
    rel_diffs = self.tester.checkJacobian("InletP0T0BC", bc_params={"phase": "0", "p0": "0.2", "T0": "0.6", "strongly_enforce_energy": True})
    n_i, n_j = rel_diffs.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)
