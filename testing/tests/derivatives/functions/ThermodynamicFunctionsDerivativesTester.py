import unittest

from ....src.testers.FunctionDerivativesTester import FunctionDerivativesTester
from sem_python.closures.thermodynamic_functions import computeSpecificVolume, computeDensity, \
  computeVelocity, computeSpecificTotalEnergy, computeSpecificInternalEnergy, \
  addKineticEnergy, subtractKineticEnergy

class ThermodynamicFunctionsFunctionDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivative_tester = FunctionDerivativesTester()

  def testVelocity(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeVelocity, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testDensity(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeDensity, 3)
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

  def testSpecificInternalEnergy(self):
    reldiffs = self.derivative_tester.checkDerivatives(computeSpecificInternalEnergy, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testAddKineticEnergy(self):
    reldiffs = self.derivative_tester.checkDerivatives(addKineticEnergy, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)

  def testSubtractKineticEnergy(self):
    reldiffs = self.derivative_tester.checkDerivatives(subtractKineticEnergy, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-7)
