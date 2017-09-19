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
from collections import OrderedDict
import warnings

# base
from enums import ModelType
from Factory import Factory

# closures
from thermodynamic_functions import computeDensity, computeVelocity, computeSpecificVolume,\
  computeSpecificTotalEnergy, computeSpecificInternalEnergy, computeSpecificEnthalpy

# input
from InputFileParser import InputFileParser
from InputFileModifier import InputFileModifier

# utilities
from error_utilities import error

## Runs the code with the given input file
# @param input_file  input file to run
# @param input_file_modifier  input file modifier object
def run(input_file, input_file_modifier=InputFileModifier()):
  # parse the input file
  input_file_parser_original = InputFileParser()
  input_file_parser_original.parse(input_file)

  # check if the input file was a differential input file
  if input_file_parser_original.blockExists("BaseInputFile"):
    # get the name of the base input file and parse it
    base_input_file = input_file_parser_original.getBlockData("BaseInputFile")["base"]
    input_file_parser_base = InputFileParser()
    input_file_parser_base.parse(base_input_file)

    # apply modifications to base input file parser
    input_file_parser_base.applyDifferentialInputFileParser(input_file_parser_original)

    input_file_parser = input_file_parser_base
  else:
    input_file_parser = input_file_parser_original

  # apply modifications to input parameters, if any
  input_file_parser.applyModifications(input_file_modifier)

  # create the factory
  factory = Factory()

  # model
  model_param_data = input_file_parser.getBlockData("Model")
  model = factory.createObject("Model", model_param_data)
  model_type = model.model_type

  # mesh(es)
  mesh_subblocks = input_file_parser.getSubblockNames("Mesh")
  if len(mesh_subblocks) == 0:
    # single mesh; no subblock is required
    mesh_param_data = input_file_parser.getBlockData("Mesh")
    mesh_class = mesh_param_data["type"]
    mesh = factory.createObject(mesh_class, mesh_param_data)
    meshes = [mesh]
  else:
    # multiple meshes; subblocks are required
    meshes = list()
    for mesh_subblock in mesh_subblocks:
      mesh_param_data = input_file_parser.getSubblockData("Mesh", mesh_subblock)
      mesh_class = mesh_param_data["type"]
      mesh_param_data["name"] = mesh_subblock
      mesh = factory.createObject(mesh_class, mesh_param_data)
      meshes.append(mesh)
  # check that no meshes have the same name
  mesh_names = list()
  for mesh in meshes:
    if mesh.name in mesh_names:
      error("Multiple meshes with the name '" + mesh.name + "'.")
    else:
      mesh_names.append(mesh.name)

  # equations of state
  eos_subblocks = input_file_parser.getSubblockNames("EoS")
  if model_type == ModelType.OnePhase:
    if len(eos_subblocks) != 1:
      error("Single-phase flow should have exactly 1 equation of state.")
  else:
    if len(eos_subblocks) != 2:
      error("Two-phase flow should have exactly 2 equations of state.")
  eos_list = list()
  phase_name_to_index = dict()
  for k, eos_subblock in enumerate(eos_subblocks):
    phase_name_to_index[eos_subblock] = k
    eos_param_data = input_file_parser.getSubblockData("EoS", eos_subblock)
    eos_class = eos_param_data["type"]
    eos_list.append(factory.createObject(eos_class, eos_param_data))

  # initial conditions / initial guess
  ic_subblocks = input_file_parser.getSubblockNames("IC")
  ics = list()
  if len(ic_subblocks) == 0:
    # ICs are to be used for every mesh
    for mesh_name in mesh_names:
      ic_param_data = input_file_parser.getBlockData("IC")
      ic_param_data["mesh_name"] = mesh_name
      if model_type == ModelType.OnePhase:
        ics.append(factory.createObject("InitialConditions1Phase", ic_param_data))
      else:
        ics.append(factory.createObject("InitialConditions2Phase", ic_param_data))
  else:
    for ic_subblock in ic_subblocks:
      ic_param_data = input_file_parser.getSubblockData("IC", ic_subblock)

      # for now, assume that IC blocks are named by corresponding mesh names
      ic_param_data["mesh_name"] = ic_subblock
      if model_type == ModelType.OnePhase:
        ics.append(factory.createObject("InitialConditions1Phase", ic_param_data))
      else:
        ics.append(factory.createObject("InitialConditions2Phase", ic_param_data))

  # DoF handler
  dof_handler_params = {"meshes": meshes, "ics": ics}
  if model_type == ModelType.OnePhase:
    dof_handler_class = "DoFHandler1Phase"
  elif model_type == ModelType.TwoPhaseNonInteracting:
    dof_handler_class = "DoFHandler2PhaseNonInteracting"
  elif model_type == ModelType.TwoPhase:
    dof_handler_class = "DoFHandler2Phase"
  dof_handler = factory.createObject(dof_handler_class, dof_handler_params)

  # junctions
  junctions = list()
  if input_file_parser.blockExists("Junctions"):
    junction_subblocks = input_file_parser.getSubblockNames("Junctions")
    for junction_subblock in junction_subblocks:
      junction_param_data = input_file_parser.getSubblockData("Junctions", junction_subblock)
      junction_class = junction_param_data["type"]
      if "phase" in junction_param_data:
        junction_param_data["phase"] = phase_name_to_index[junction_param_data["phase"]]
      junction_param_data["dof_handler"] = dof_handler
      junction_param_data["eos_list"] = eos_list

      junction = factory.createObject(junction_class, junction_param_data)
      junctions.append(junction)

  # update DoF handler with junction constraints
  dof_handler.updateWithJunctionConstraints(junctions)

  # boundary conditions
  bcs = list()
  bc_subblocks = input_file_parser.getSubblockNames("BC")
  for bc_subblock in bc_subblocks:
    bc_param_data = input_file_parser.getSubblockData("BC", bc_subblock)
    bc_class = bc_param_data["type"]

    # if there is only 1 mesh, add the mesh parameter if not supplied already
    if len(meshes) == 1:
      if "mesh_name" not in bc_param_data:
        bc_param_data["mesh_name"] = meshes[0].name

    if "phase" in bc_param_data:
      bc_param_data["phase"] = phase_name_to_index[bc_param_data["phase"]]
    bc_param_data["dof_handler"] = dof_handler
    bc_param_data["eos_list"] = eos_list

    bc = factory.createObject(bc_class, bc_param_data)
    bcs.append(bc)

  # interface closures
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
  if len(gravity) != 3:
    error("Gravity vector must have 3 elements")

  # heat transfer data
  # initialize heat transfer data such that no heat transfer occurs
  ht_param_data_default = {"T_wall": 0, "htc_wall": 0, "P_heat": 1}
  ht_data_default = factory.createObject("HeatTransferData", ht_param_data_default)
  ht_data = [ht_data_default] * len(meshes)
  ht_subblocks = input_file_parser.getSubblockNames("HeatTransfer")
  for subblock in ht_subblocks:
    # sub-blocks should be named by the corresponding mesh names
    mesh_name = subblock
    i_mesh = dof_handler.mesh_name_to_mesh_index[mesh_name]
    ht_data_mesh = input_file_parser.getSubblockData("HeatTransfer", subblock)
    ht_data[i_mesh] = factory.createObject("HeatTransferData", ht_data_mesh)

  # nonlinear solver options
  nonlinear_solver_params = input_file_parser.getBlockData("NonlinearSolver")

  # stabilization
  if input_file_parser.blockExists("Stabilization"):
    stabilization_param_data = input_file_parser.getBlockData("Stabilization")
    stabilization_param_data["factory"] = factory
    stabilization_param_data["dof_handler"] = dof_handler
    stabilization_param_data["model_type"] = model_type
    stabilization_class = stabilization_param_data["type"]
  else:
    stabilization_param_data = {"factory": factory, "dof_handler": dof_handler, "model_type": model_type}
    stabilization_class = "NoStabilization"
  stabilization = factory.createObject(stabilization_class, stabilization_param_data)

  # create and run the executioner
  executioner_param_data = input_file_parser.getBlockData("Executioner")
  executioner_type = executioner_param_data["type"]
  executioner_param_data["model"] = model
  executioner_param_data["ics"] = ics
  executioner_param_data["bcs"] = bcs
  executioner_param_data["junctions"] = junctions
  executioner_param_data["eos_list"] = eos_list
  executioner_param_data["interface_closures"] = interface_closures
  executioner_param_data["gravity"] = gravity
  executioner_param_data["dof_handler"] = dof_handler
  executioner_param_data["meshes"] = meshes
  executioner_param_data["nonlinear_solver_params"] = nonlinear_solver_params
  executioner_param_data["stabilization"] = stabilization
  executioner_param_data["factory"] = factory
  executioner = factory.createObject(executioner_type, executioner_param_data)
  U = executioner.run()

  # perform post-processing

  if input_file_parser.blockExists("Output"):
    # get solution and compute aux quantities
    vf1, arhoA1, arhouA1, arhoEA1 = dof_handler.getPhaseSolution(U, 0)
    if (model_type != ModelType.OnePhase):
      vf2, arhoA2, arhouA2, arhoEA2 = dof_handler.getPhaseSolution(U, 1)
    rho1 = computeDensity(vf1, arhoA1, dof_handler.A)[0]
    u1 = computeVelocity(arhoA1, arhouA1)[0]
    v1 = computeSpecificVolume(rho1)[0]
    E1 = computeSpecificTotalEnergy(arhoA1, arhoEA1)[0]
    e1 = computeSpecificInternalEnergy(u1, E1)[0]
    eos1 = eos_list[0]
    p1 = eos1.p(v1, e1)[0]
    T1 = eos1.T(v1, e1)[0]
    h1 = computeSpecificEnthalpy(e1, p1, rho1)[0]
    H1 = h1 + 0.5 * u1**2
    s1 = eos1.s(v1, e1)[0]
    # overflow can occur in p(h,s) with small gamma values
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      p01 = eos1.p_from_h_s(H1, s1)[0]
    if (model_type != ModelType.OnePhase):
      rho2 = computeDensity(vf2, arhoA2, dof_handler.A)[0]
      u2 = computeVelocity(arhoA2, arhouA2)[0]
      v2 = computeSpecificVolume(rho2)[0]
      E2 = computeSpecificTotalEnergy(arhoA2, arhoEA2)[0]
      e2 = computeSpecificInternalEnergy(u2, E2)[0]
      eos2 = eos_list[1]
      p2 = eos2.p(v2, e2)[0]
      T2 = eos2.T(v2, e2)[0]
      h2 = computeSpecificEnthalpy(e2, p2, rho2)[0]
      H2 = h2 + 0.5 * u2**2
      s2 = eos2.s(v2, e2)[0]
      # overflow can occur in p(h,s) with small gamma values
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p02 = eos2.p_from_h_s(H2, s2)[0]

    # create an ordered data dictionary
    data = OrderedDict()
    data["x"] = dof_handler.x
    data["A"] = dof_handler.A
    if model_type == ModelType.OnePhase:
      data["rhoA"] = arhoA1
      data["rhouA"] = arhouA1
      data["rhoEA"] = arhoEA1
      data["rho"] = rho1
      data["u"] = u1
      data["v"] = v1
      data["E"] = E1
      data["e"] = e1
      data["p"] = p1
      data["T"] = T1
      data["h"] = h1
      data["H"] = H1
      data["s"] = s1
      data["p0"] = p01
    else:
      # phase 1
      data["vf1"] = vf1
      data["arhoA1"] = arhoA1
      data["arhouA1"] = arhouA1
      data["arhoEA1"] = arhoEA1
      data["rho1"] = rho1
      data["u1"] = u1
      data["v1"] = v1
      data["E1"] = E1
      data["e1"] = e1
      data["p1"] = p1
      data["T1"] = T1
      data["h1"] = h1
      data["H1"] = H1
      data["s1"] = s1
      data["p01"] = p01

      # phase 2
      data["vf2"] = vf2
      data["arhoA2"] = arhoA2
      data["arhouA2"] = arhouA2
      data["arhoEA2"] = arhoEA2
      data["rho2"] = rho2
      data["u2"] = u2
      data["v2"] = v2
      data["E2"] = E2
      data["e2"] = e2
      data["p2"] = p2
      data["T2"] = T2
      data["h2"] = h2
      data["H2"] = H2
      data["s2"] = s2
      data["p02"] = p02

    if model_type == ModelType.OnePhase:
      default_print_list = ["rho", "u", "p"]
    elif model_type == ModelType.TwoPhaseNonInteracting:
      default_print_list = ["rho1", "u1", "p1", "rho2", "u2", "p2"]
    else: # model_type == ModelType.TwoPhase
      default_print_list = ["vf1", "rho1", "u1", "p1", "rho2", "u2", "p2"]

    # create and run outputs
    output_subblocks = input_file_parser.getSubblockNames("Output")
    for output_subblock in output_subblocks:
      # extract the output parameter data
      output_param_data = input_file_parser.getSubblockData("Output", output_subblock)
      output_class = output_param_data["type"]
      output_param_data["dof_handler"] = dof_handler
      # add additional parameters for specific classes
      if output_class == "ScreenOutput":
        if "data_names" not in output_param_data:
          output_param_data["data_names"] = default_print_list
      # create the output
      output = factory.createObject(output_class, output_param_data)

      # run the output
      output.run(data)

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
