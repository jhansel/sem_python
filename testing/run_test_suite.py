# to get test results in color, install "colour-runner":
#   pip install colour-runner
try:
  from colour_runner.runner import ColourTextTestRunner
  print_in_color = True
except:
  from unittest import TextTestRunner
  print_in_color = False

import os
import re
import sys
from unittest import TestSuite, defaultTestLoader

# add source file paths
source_directories = os.listdir("../src")
for source_dir in source_directories:
  sys.path.append("../src/" + source_dir)

# tests
for root, subdirs, files in os.walk("tests"):
  for subdir in subdirs:
    sys.path.append("/".join([root, subdir]))

# testing code
sys.path.append("src/utilities")

def getTestModuleList():
  module_regex = re.compile(r'.*\.py')

  test_modules = list()
  for root, subdirs, files in os.walk("tests"):
    for thefile in files:
      fields = thefile.split(".")
      ext = fields[-1]
      if ext == "py":
        test_modules.append(".".join(fields[0:-1]))

  return test_modules

def main():
  # list of test modules
  test_modules = getTestModuleList()

  # add tests to test suite
  suite = TestSuite()
  for test_module in test_modules:
    suite.addTest(defaultTestLoader.loadTestsFromName(test_module))

  # run test suite
  if print_in_color:
    ColourTextTestRunner(verbosity=2).run(suite)
  else:
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
  main()
