import sys
sys.path.append("../../src/utilities")
sys.path.append("../../../src/base")
sys.path.append("../../../src/closures")
sys.path.append("../../../src/input")
sys.path.append("../../../src/utilities")

import unittest

from DerivativesTester import DerivativesTester
from IdealGasEoS import IdealGasEoS, IdealGasEoSParameters

class IdealGasEoSDerivativesTester(unittest.TestCase):
  def setUp(self):
    params = IdealGasEoSParameters()
    params.set("gamma", 1.4)
    params.set("R", 290.0)
    self.eos = IdealGasEoS(params)
    self.derivative_tester = DerivativesTester()

  def testSpecificInternalEnergy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.e, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testPressure(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testTemperature(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.T, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

if __name__ == "__main__":
  params = IdealGasEoSParameters()
  params.set("gamma", 1.4)
  params.set("R", 290.0)

  eos = IdealGasEoS(params)

  tester = DerivativesTester(False)
  tester.checkDerivatives(eos.e, 2)
  tester.checkDerivatives(eos.p, 2)
  tester.checkDerivatives(eos.T, 2)
