import unittest

from FunctionDerivativesTester import FunctionDerivativesTester
from StiffenedGasEoS import StiffenedGasEoS, StiffenedGasEoSParameters

class StiffenedGasEoSDerivativesTester(unittest.TestCase):
  def setUp(self):
    params = StiffenedGasEoSParameters()
    params.set("gamma", 1.4)
    params.set("cv", 2000)
    params.set("q", -5)
    params.set("p_inf", 10)
    self.eos = StiffenedGasEoS(params)
    self.derivative_tester = FunctionDerivativesTester()

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

  def testSoundSpeed(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.c, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

if __name__ == "__main__":
  params = StiffenedGasEoSParameters()
  params.set("gamma", 1.4)
  params.set("cv", 2000)
  params.set("q", -5)
  params.set("p_inf", 10)

  eos = StiffenedGasEoS(params)

  tester = FunctionDerivativesTester(False)
  tester.checkDerivatives(eos.e, 2)
  tester.checkDerivatives(eos.p, 2)
  tester.checkDerivatives(eos.T, 2)
  tester.checkDerivatives(eos.c, 2)
