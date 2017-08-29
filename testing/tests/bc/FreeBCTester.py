import unittest

from BCTester import BCTester

class FreeBCTester(unittest.TestCase):
  def setUp(self):
    self.tester = BCTester()

  def testJacobian(self):
    rel_diffs = self.tester.checkJacobian("FreeBC")
    for var_i_dict in rel_diffs:
      for var_j in var_i_dict:
        self.assertLessEqual(var_i_dict[var_j], 1e-6)

if __name__ == "__main__":
  tester = BCTester(True)
  _ = tester.checkJacobian("FreeBC")
