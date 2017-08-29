from copy import deepcopy
import numpy as np

from display_utilities import computeRelativeDifferenceMatrix, printRelativeMatrixDifference
from enums import ModelType
from error_utilities import error
from Factory import Factory
from numeric_utilities import computeRelativeDifference

class JunctionTester(object):
  def __init__(self, junction_name, verbose=False):
    self.junction_name = junction_name
    self.verbose = verbose

  def checkJacobian(self, test_option, model_type=ModelType.OnePhase, phase=0, junction_params=dict(), fd_eps=1e-8):
    # factory
    factory = Factory()

    # meshes
    params1 = {"n_cell": 1, "name": "mesh1"}
    params2 = {"n_cell": 1, "name": "mesh2"}
    meshes = [factory.createObject("UniformMesh", params1), factory.createObject("UniformMesh", params2)]

    # DoF handler
    dof_handler_params = {"meshes": meshes}
    if model_type == ModelType.OnePhase:
      dof_handler_class = "DoFHandler1Phase"
    elif model_type == ModelType.TwoPhaseNonInteracting:
      dof_handler_class = "DoFHandler2PhaseNonInteracting"
      def vf1_initial(x):
        return 0.3
      dof_handler_params["initial_vf1"] = vf1_initial
    elif model_type == ModelType.TwoPhase:
      dof_handler_class = "DoFHandler2Phase"
    dof_handler = factory.createObject(dof_handler_class, dof_handler_params)
    n_dof = dof_handler.n_dof
    n_var = dof_handler.n_var

    # equation of state
    eos_params1 = {"slope_initial": 1.0, "slope_increment": 0.1}
    eos_list = [factory.createObject("TestEoS", eos_params1)]
    if model_type != ModelType.OnePhase:
      eos_params2 = {"slope_initial": 1.1, "slope_increment": 0.07}
      eos_list.append(factory.createObject("TestEoS", eos_params2))

    # junction
    junction_params["mesh_names"] = " ".join([mesh.name for mesh in meshes])
    junction_params["mesh_sides"] = "right left"
    junction_params["dof_handler"] = dof_handler
    junction_params["eos_list"] = eos_list
    junction_parameters = factory.createParametersObject(self.junction_name)
    for param in junction_params:
      junction_parameters.set(param, junction_params[param])
    if junction_parameters.hasRegisteredParam("phase"):
      junction_parameters.set("phase", phase)
    junction = factory.createObjectFromParametersObject(self.junction_name, junction_parameters)

    # compute base solution
    U = np.zeros(n_dof)
    U_old = np.zeros(n_dof)
    for i in xrange(n_dof):
      U[i] = i + 1.0
      U_old[i] = i + 2.0

    # determine evaluation function
    if test_option == "weak":
      f = junction.applyWeaklyToNonlinearSystem
    elif test_option == "strong":
      f = junction.applyStronglyToNonlinearSystem
    elif test_option == "both":
      def f(*args, **kwargs):
        junction.applyWeaklyToNonlinearSystem(*args, **kwargs)
        junction.applyStronglyToNonlinearSystem(*args, **kwargs)
    else:
      error("Invalid test option")

    # base calculation
    r = np.zeros(n_dof)
    J_hand_coded = np.zeros(shape=(n_dof, n_dof))
    f(U, U_old, r, J_hand_coded)

    # finite difference Jacobians
    rel_diffs = np.zeros(shape=(n_dof, n_dof))
    J_fd = np.zeros(shape=(n_dof, n_dof))
    for j in xrange(n_dof):
      # perturb solution
      U_perturbed = deepcopy(U)
      U_perturbed[j] += fd_eps

      # compute finite difference Jacobian
      r_perturbed = np.zeros(n_dof)
      J_perturbed = np.zeros(shape=(n_dof, n_dof))
      f(U_perturbed, U_old, r_perturbed, J_perturbed)
      for i in xrange(n_dof):
        J_fd[i,j] = (r_perturbed[i] - r[i]) / fd_eps

    # compute relative difference matrix
    rel_diffs = computeRelativeDifferenceMatrix(J_hand_coded, J_fd)

    # print results
    if self.verbose:
      print "\nRelative difference of Jacobian for " + test_option + " contributions:"
      printRelativeMatrixDifference(rel_diffs, J_hand_coded - J_fd, 1e-1, 1e-3)

    # take the absolute value of the relative differences
    return abs(rel_diffs)
