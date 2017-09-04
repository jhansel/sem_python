import unittest
from numpy import sin, cos, exp, log

from sem_python.utilities.conversion_utilities import stringToFunction

class StringToFunctionTester(unittest.TestCase):
  def setUp(self):
    f_string = "sin(exp(x)) if x < 1 else cos(log(x))"
    self.f = stringToFunction(f_string)

  def testIf(self):
    y = 0.5
    result_by_hand = sin(exp(y))
    result_from_f = self.f(y)
    self.assertEqual(result_by_hand, result_from_f)

  def testElse(self):
    y = 1.5
    result_by_hand = cos(log(y))
    result_from_f = self.f(y)
    self.assertEqual(result_by_hand, result_from_f)
