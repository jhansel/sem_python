import unittest

from FunctionDerivativesTester import FunctionDerivativesTester
from StiffenedGasEoS import StiffenedGasEoS, StiffenedGasEoSParameters

class StiffenedGasEoSDerivativesTester(unittest.TestCase):
  def setUp(self):
    params = StiffenedGasEoSParameters()
    params.set("gamma", 1.4)
    params.set("cv", 2.0)
    params.set("q", -5.0)
    params.set("p_inf", 0.75)
    params.set("q_prime", 6.0)
    self.eos = StiffenedGasEoS(params)
    self.derivative_tester = FunctionDerivativesTester()

  def testDensity(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.rho, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testDensityFromPressureEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.rho, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

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
      self.assertLessEqual(reldiff, 5e-6)

  def testPressureFromEnthalpyEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p_from_h_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

if __name__ == "__main__":
  params = StiffenedGasEoSParameters()
  params.set("gamma", 1.4)
  params.set("cv", 2.0)
  params.set("q", -5.0)
  params.set("p_inf", 0.75)
  params.set("q_prime", 6.0)

  eos = StiffenedGasEoS(params)

  tester = FunctionDerivativesTester(False)
  tester.checkDerivatives(eos.rho, 2)
  tester.checkDerivatives(eos.rho_from_p_s, 2)
  tester.checkDerivatives(eos.e, 2)
  tester.checkDerivatives(eos.p, 2)
  tester.checkDerivatives(eos.T, 2)
  tester.checkDerivatives(eos.c, 2)
  tester.checkDerivatives(eos.s, 2)
  tester.checkDerivatives(eos.p_from_h_s, 2)
