import unittest

import sem
from ....src.testers.JunctionTester import JunctionTester

class EqualFluxLM1PhaseJunctionTester(unittest.TestCase):
  # Currently the LM contributions to the non-constraint equations still need
  # to be derived and implemented, so this test is disabled for now.
  def testDummy(self):
    self.assertTrue(True)

if __name__ == "__main__":
  tester = JunctionTester("EqualFluxLM1PhaseJunction", verbose=True)
  _ = tester.checkJacobian("both")
