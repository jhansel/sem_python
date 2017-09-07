import unittest

from FunctionDerivativesTester import FunctionDerivativesTester
from IdealGasEoS import IdealGasEoS, IdealGasEoSParameters

class IdealGasEoSFunctionDerivativesTester(unittest.TestCase):
  def setUp(self):
    params = IdealGasEoSParameters()
    params.set("gamma", 1.4)
    params.set("R", 290.0)
    self.eos = IdealGasEoS(params)
    self.derivative_tester = FunctionDerivativesTester()

  def testDensity(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.rho, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testDensityFromPressureEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.rho_from_p_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 5e-6)

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

  def testEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testPressureFromEnthalpyEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p_from_h_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 5e-6)

if __name__ == "__main__":
  params = IdealGasEoSParameters()
  params.set("gamma", 1.4)
  params.set("R", 290.0)

  eos = IdealGasEoS(params)

  tester = FunctionDerivativesTester(False)
  tester.checkDerivatives(eos.e, 2)
  tester.checkDerivatives(eos.p, 2)
  tester.checkDerivatives(eos.T, 2)
  tester.checkDerivatives(eos.c, 2)
  tester.checkDerivatives(eos.s, 2)
  tester.checkDerivatives(eos.p_from_h_s, 2)
