import unittest

import sem
from sem_python.input.InputFileModifier import InputFileModifier
from ....src.testers.CSVTester import CSVTester
from ....src.testers.JunctionTester import JunctionTester

class EqualFluxJunctionTester(unittest.TestCase):
  def testJacobian(self):
    tester = JunctionTester("EqualFluxJunction")
    matched = tester.checkJacobian("both")
    n_i, n_j = matched.shape
    for i in xrange(n_i):
      for j in xrange(n_j):
        self.assertTrue(matched[i,j])
