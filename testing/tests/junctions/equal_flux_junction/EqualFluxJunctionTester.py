import unittest

import sem
from sem_python.input.ParameterModification import BlockParameterModification
from .. import CSVTester
from ..JunctionTester import JunctionTester

class EqualFluxJunctionTester(unittest.TestCase):
  def testJacobian(self):
    tester = JunctionTester("EqualFluxJunction")
    rel_diffs = tester.checkJacobian("both")
    n_i, n_j = rel_diffs.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertLessEqual(rel_diffs[i,j], 1e-6)

if __name__ == "__main__":
  tester = JunctionTester("EqualFluxJunction", verbose=True)
  _ = tester.checkJacobian("both")
