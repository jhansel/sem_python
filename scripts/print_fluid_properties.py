"""
This script prints a number of fluid properties at a given state.

Usage:
  print_fluid_properties.py -h | --help

Options:
  -h --help  Display help information.
"""

from docopt import docopt

from sem_python.base.Factory import Factory
from sem_python.closures.ThermodynamicState import ThermodynamicState
from sem_python.input.InputFileParser import InputFileParser

## Runs the script
# @param mods  input file modifications
def run(mods=list()):
  # parse the input file
  input_file_parser = InputFileParser()
  input_file_parser.parse("print_fluid_properties.in")

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

if __name__ == "__main__":
  run()
