import unittest

from ...src.testers.FunctionDerivativesTester import FunctionDerivativesTester
from sem_python.closures.TestEoS import TestEoS, TestEoSParameters

class TestEoSTester(unittest.TestCase):
  def setUp(self):
    params = TestEoSParameters()
    self.eos = TestEoS(params)
    self.derivative_tester = FunctionDerivativesTester()

  def testSpecificInternalEnergyDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.e, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testPressureDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testTemperatureDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.T, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testSoundSpeedDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.c, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testEntropyDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testEntropyFromEnthalpyPressureDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.s_from_h_p, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 5e-6)

  def testPressureFromEnthalpyEntropyDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p_from_h_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testSpecificEnthalpyDerivatives(self):
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
  tester.checkDerivatives(eos.s_from_h_p, 2)
  tester.checkDerivatives(eos.p_from_h_s, 2)
