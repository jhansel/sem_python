import sys

from file_utilities import readCSVData

class CSVTester(object):
  def __init__(self, test_dir, test_file, abs_tol=1e-12):
    self.test_dir = test_dir
    self.test_file = test_file
    self.abs_tol = abs_tol

  def filesAreEqual(self):
    # create full file paths
    if self.test_dir.endswith("/"):
      test_file = self.test_dir + self.test_file
      gold_file = self.test_dir + "gold/" + self.test_file
    else:
      test_file = self.test_dir + "/" + self.test_file
      gold_file = self.test_dir + "/gold/" + self.test_file

    # read the data from both files
    test_data = readCSVData(test_file)
    gold_data = readCSVData(gold_file)

    # exit if keys are not equal
    if test_data.keys() != gold_data.keys():
      print "\nDifferent keys found:"
      print "Gold keys: ", gold_data.keys()
      print "Test keys: ", test_data.keys()
      return False

    # loop over keys
    files_are_equal = True
    first_diff_encountered = False
    for key in test_data:
      y_test = test_data[key]
      y_gold = gold_data[key]

      # exit if sizes are not equal
      if len(y_test) != len(y_gold):
        print "\nDifferent sizes: gold size = " + len(y_gold) + ", test size = " + len(y_test)
        return False

      # loop over all values
      for i, value in enumerate(y_test):
        if abs(value - y_gold[i]) > self.abs_tol:
          if not first_diff_encountered:
            print "" # print newline for formatting in test suite
            first_diff_encountered = True
          print "Diff: %s[%i]: gold = %13.6e, test = %13.6e, diff = %13.6e" % (key, i, y_gold[i], value, y_gold[i] - value)
          files_are_equal = False

    return files_are_equal
