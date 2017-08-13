"""
This script solves the seven-equation model (SEM) of 1-D two-phase flow
using the Continuous finite element method (CFEM). For help, please see
the external documentation in doc/.

Usage:
  sem.py -h | --help
  sem.py <INPUTFILE>

Options:
  -h --help  Display help information.

Arguments:
  INPUTFILE  Name of the input file
"""

from docopt import docopt

# base
from enums import ModelType
from Factory import Factory

# input
from InputFileParser import InputFileParser

# utilities
from error_utilities import error

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

  # model
  model_param_data = input_file_parser.getBlockData("Model")
  model = factory.createObject("Model", model_param_data)
  model_type = model.model_type

  # mesh
  mesh_param_data = input_file_parser.getBlockData("Mesh")
  mesh_class = mesh_param_data["type"]
  mesh = factory.createObject(mesh_class, mesh_param_data)

  # equations of state
  eos_subblocks = input_file_parser.getSubblockNames("EoS")
  if model_type == ModelType.OnePhase:
    if len(eos_subblocks) != 1:
      error("Model type '" + model_type + "' should have exactly 1 equation of state.")
  else:
    if len(eos_subblocks) != 2:
      error("Model type '" + model_type + "' should have exactly 2 equations of state.")
  eos = list()
  phase_name_to_index = dict()
  for k, eos_subblock in enumerate(eos_subblocks):
    phase_name_to_index[eos_subblock] = k
    eos_param_data = input_file_parser.getSubblockData("EoS", eos_subblock)
    eos_class = eos_param_data["type"]
    eos.append(factory.createObject(eos_class, eos_param_data))

  # initial conditions / initial guess
  ic_param_data = input_file_parser.getBlockData("IC")
  if model_type == ModelType.OnePhase:
    ics = factory.createObject("InitialConditions1Phase", ic_param_data)
  else:
    ics = factory.createObject("InitialConditions2Phase", ic_param_data)

  # DoF handler
  dof_handler_params = {"mesh": mesh}
  if model_type == ModelType.OnePhase:
    dof_handler_class = "DoFHandler1Phase"
  elif model_type == ModelType.TwoPhaseNonInteracting:
    dof_handler_class = "DoFHandler2PhaseNonInteracting"
    dof_handler_params["initial_vf1"] = ics.vf1
  elif model_type == ModelType.TwoPhase:
    dof_handler_class = "DoFHandler2Phase"
  dof_handler = factory.createObject(dof_handler_class, dof_handler_params)

  # boundary conditions
  bcs = list()
  bc_subblocks = input_file_parser.getSubblockNames("BC")
  for bc_subblock in bc_subblocks:
    bc_param_data = input_file_parser.getSubblockData("BC", bc_subblock)
    bc_class = bc_param_data["type"]
    if "phase" in bc_param_data:
      bc_param_data["phase"] = phase_name_to_index[bc_param_data["phase"]]
    bc_param_data["dof_handler"] = dof_handler
    bc_param_data["eos"] = eos
    bc = factory.createObject(bc_class, bc_param_data)
    bcs.append(bc)

  # interace closures
  if model_type == ModelType.TwoPhase:
    interface_closures_params = input_file_parser.getBlockData("InterfaceClosures")
    interface_closures_params["factory"] = factory
    interface_closures_class = interface_closures_params["type"]
    interface_closures = factory.createObject(interface_closures_class, interface_closures_params)
  else:
    interface_closures = None

  # gravity
  physics_param_data = input_file_parser.getBlockData("Physics")
  physics_params = factory.createParametersObject("Physics", physics_param_data)
  gravity = physics_params.get("gravity")

  # nonlinear solver options
  nonlinear_solver_params = input_file_parser.getBlockData("NonlinearSolver")

  # stabilization
  if input_file_parser.blockExists("Stabilization"):
    stabilization_param_data = input_file_parser.getBlockData("Stabilization")
    stabilization_param_data["factory"] = factory
    stabilization_param_data["dof_handler"] = dof_handler
    stabilization_class = stabilization_param_data["type"]
  else:
    stabilization_param_data = {"factory": factory, "dof_handler": dof_handler}
    stabilization_class = "NoStabilization"
  stabilization = factory.createObject(stabilization_class, stabilization_param_data)

  # create and run the executioner
  executioner_param_data = input_file_parser.getBlockData("Executioner")
  executioner_type = executioner_param_data["type"]
  executioner_param_data["model"] = model
  executioner_param_data["ics"] = ics
  executioner_param_data["bcs"] = bcs
  executioner_param_data["eos"] = eos
  executioner_param_data["interface_closures"] = interface_closures
  executioner_param_data["gravity"] = gravity
  executioner_param_data["dof_handler"] = dof_handler
  executioner_param_data["mesh"] = mesh
  executioner_param_data["nonlinear_solver_params"] = nonlinear_solver_params
  executioner_param_data["stabilization"] = stabilization
  executioner_param_data["factory"] = factory
  executioner = factory.createObject(executioner_type, executioner_param_data)
  U = executioner.run()

  # create and run the postprocessor
  postprocessor_param_data = input_file_parser.getBlockData("Output")
  postprocessor_param_data["model"] = model
  postprocessor_param_data["eos"] = eos
  postprocessor_param_data["dof_handler"] = dof_handler
  postprocessor_param_data["mesh"] = mesh
  postprocessor = factory.createObject("Postprocessor", postprocessor_param_data)
  postprocessor.run(U)

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
