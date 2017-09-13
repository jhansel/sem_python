import unittest

from FunctionDerivativesTester import FunctionDerivativesTester
from TestEoS import TestEoS, TestEoSParameters

class TestEoSFunctionDerivativesTester(unittest.TestCase):
  def setUp(self):
    params = TestEoSParameters()
    self.eos = TestEoS(params)
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

  def testEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testPressureFromEnthalpyEntropy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p_from_h_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testSpecificEnthalpy(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.h, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

if __name__ == "__main__":
  params = TestEoSParameters()
  eos = TestEoS(params)

  tester = FunctionDerivativesTester(False)
  tester.checkDerivatives(eos.e, 2)
  tester.checkDerivatives(eos.p, 2)
  tester.checkDerivatives(eos.T, 2)
  tester.checkDerivatives(eos.c, 2)
  tester.checkDerivatives(eos.s, 2)
  tester.checkDerivatives(eos.p_from_h_s, 2)
