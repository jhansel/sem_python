import unittest

from FunctionDerivativesTester import FunctionDerivativesTester
from thermodynamic_functions import computeSpecificVolume, computeDensity, \
  computeVelocity, computeSpecificTotalEnergy, computeSpecificInternalEnergy

class ThermodynamicFunctionsFunctionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivative_tester = FunctionDerivativesTester()

  def testVelocity(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeVelocity, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testDensity(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeDensity, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testSpecificVolume(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeSpecificVolume, 1)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testSpecificTotalEnergy(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeSpecificTotalEnergy, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testSpecificInternalEnergy(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeSpecificInternalEnergy, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

if __name__ == "__main__":
  tester = FunctionDerivativesTester(False)
  tester.checkDerivatives(computeVelocity, 2)
  tester.checkDerivatives(computeDensity, 2)
  tester.checkDerivatives(computeSpecificVolume, 1)
  tester.checkDerivatives(computeSpecificTotalEnergy, 2)
  tester.checkDerivatives(computeSpecificInternalEnergy, 2)
