"""
This script prints a number of fluid properties at a given state.

Usage:
  print_fluid_properties.py -h | --help
  print_fluid_properties.py <INPUTFILE>

Options:
  -h --help  Display help information.

Arguments:
  INPUTFILE  Name of input file
"""

from docopt import docopt

from Factory import Factory
from ThermodynamicState import ThermodynamicState
from InputFileParser import InputFileParser

## Runs the code with the given input file
# @param input_file  input file to run
# @param mods  input file modifications
def run(input_file, mods=list()):
  # parse the input file
  input_file_parser = InputFileParser()
  input_file_parser.parse(input_file)

  # apply modifications to input parameters, if any
  for mod in mods:
    input_file_parser.applyModification(mod)

  # create the factory
  factory = Factory()

  # equation of state
  eos_param_data = input_file_parser.getBlockData("EoS")
  eos_class = eos_param_data["type"]
  eos = factory.createObject(eos_class, eos_param_data)

  # thermodynamic state
  state_data = input_file_parser.getBlockData("ThermodynamicState")
  state = factory.createObject("ThermodynamicState", state_data)
  state.computeRemainingProperties(eos)
  print state

## Gets the input file from the command line and runs it
def main():
  # parse command-line arguments using docopt
  command_line_arguments = docopt(__doc__)

  # get the name of the input file from the command line
  input_file = command_line_arguments["<INPUTFILE>"]

  # run with the input file
  run(input_file)

if __name__ == "__main__":
  main()
