import unittest

from ...src.testers.EoSConsistencyTester import EoSConsistencyTester
from ...src.testers.FunctionDerivativesTester import FunctionDerivativesTester
from sem_python.closures.IdealGasEoS import IdealGasEoS, IdealGasEoSParameters

class IdealGasEoSTester(unittest.TestCase):
  def setUp(self):
    params = IdealGasEoSParameters()
    params.set("gamma", 1.4)
    params.set("R", 290.0)
    self.eos = IdealGasEoS(params)
    self.derivative_tester = FunctionDerivativesTester()

  def testDensityDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.rho, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testDensityFromPressureEntropyDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.rho_from_p_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 5e-6)

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

  def testPressureFromEnthalpyEntropyDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.p_from_h_s, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 5e-6)

  def testSpecificEnthalpyDerivatives(self):
    reldiffs = self.derivative_tester.checkDerivatives(self.eos.h, 2)
    for reldiff in reldiffs:
      self.assertLessEqual(reldiff, 1e-6)

  def testEoSConsistency(self):
    eos_consistency_tester = EoSConsistencyTester(False)
    reldiffs = eos_consistency_tester.checkConsistency(self.eos)
    for check in reldiffs:
      self.assertLessEqual(reldiffs[check], 1e-12)

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

  eos_consistency_tester = EoSConsistencyTester(True)
  _ = eos_consistency_tester.checkConsistency(eos)
