import unittest

from sem_python.utilities.string_utilities import stripStringFromRight

class StringUtilitiesTester(unittest.TestCase):
  def testStripStringFromRightChanged(self):
    self.assertTrue(stripStringFromRight("mystrings.csv", ".csv") == "mystrings")
  def testStripStringFromRightUnchanged(self):
    self.assertTrue(stripStringFromRight("mystrings.pdf", ".csv") == "mystrings.pdf")
