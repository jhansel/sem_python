"""
This script runs the test suite or a subset of the test suite.

Usage:
  run_tests.py -h | --help
  run_tests.py [<REGEX>]

Options:
  -h --help  Display help information.

Arguments:
  REGEX  Regular expression for test modules to include
"""

# To get test results in color, install "colour-runner":
#   pip install colour-runner
try:
  from colour_runner.runner import ColourTextTestRunner
  print_in_color = True
except:
  from unittest import TextTestRunner
  print_in_color = False

from docopt import docopt
import os
import re
import sys
from unittest import defaultTestLoader

def main():
  # parse command-line arguments using docopt
  command_line_arguments = docopt(__doc__)

  # get the the regex parameter, if any
  regex = command_line_arguments["<REGEX>"]
  if regex:
    search_pattern = "*" + regex + "*Tester.py"
  else:
    search_pattern = "*Tester.py"

  # add tests to test suite
  suite = defaultTestLoader.discover(
    start_dir='testing/tests', pattern=search_pattern, top_level_dir='.')

  # run test suite
  if print_in_color:
    ColourTextTestRunner(verbosity=2).run(suite)
  else:
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
  main()
