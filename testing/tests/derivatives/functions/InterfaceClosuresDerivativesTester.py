import unittest

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/testing/src")
from FunctionDerivativesTester import FunctionDerivativesTester

sys.path.append(base_dir + "src/closures")
from InterfaceClosures import InterfaceClosures, InterfaceClosuresParameters

class InterfaceClosuresFunctionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivative_tester = FunctionDerivativesTester()
    params = InterfaceClosuresParameters()
    params.set("chi", 0.5)
    params.set("pressure_relaxation_time", 0.1)
    self.closures = InterfaceClosures(params)

  def testInterfaceVelocity(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.closures.computeInterfaceVelocity, 3)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testInterfacePressure(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.closures.computeInterfacePressure, 3)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testBeta(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.closures.computeBeta, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testMu(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.closures.computeMu, 3)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 5e-7)

  def testTheta(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.closures.computeTheta, 3)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

if __name__ == "__main__":
  params = InterfaceClosuresParameters()
  params.set("chi", 0.5)
  params.set("pressure_relaxation_time", 0.1)
  closures = InterfaceClosures(params)

  tester = FunctionDerivativesTester(False)
  tester.checkDerivatives(closures.computeInterfaceVelocity, 3)
  tester.checkDerivatives(closures.computeInterfacePressure, 3)
  tester.checkDerivatives(closures.computeBeta, 2)
  tester.checkDerivatives(closures.computeMu, 3)
  tester.checkDerivatives(closures.computeTheta, 3)
