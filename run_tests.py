# To get test results in color, install "colour-runner":
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
from unittest import defaultTestLoader

def main():
  # add tests to test suite
  suite = defaultTestLoader.discover(
    start_dir='testing/tests', pattern='*Tester.py', top_level_dir='.')

  # run test suite
  if print_in_color:
    ColourTextTestRunner(verbosity=2).run(suite)
  else:
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
  main()
