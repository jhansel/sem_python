from copy import deepcopy
import numpy as np

from sem_python.utilities.display_utilities import computeRelativeDifferenceMatrix, printRelativeMatrixDifference, printMatrix
from sem_python.base.enums import ModelType
from sem_python.utilities.error_utilities import error
from sem_python.base.Factory import Factory
from sem_python.utilities.numeric_utilities import computeRelativeDifference

class JunctionTester(object):
  def __init__(self, junction_name, rel_tol=1e-6, abs_tol=1e-6):
    self.junction_name = junction_name
    self.rel_tol = rel_tol
    self.abs_tol = abs_tol

  def checkJacobian(self, test_option, model_type=ModelType.OnePhase, phase=0,
      junction_params=dict(), fd_eps=1e-8, use_zero_velocity=False):
    # factory
    factory = Factory()

    # meshes
    params1 = {"n_cell": 1, "name": "mesh1"}
    params2 = {"n_cell": 1, "name": "mesh2"}
    meshes = [factory.createObject("UniformMesh", params1), factory.createObject("UniformMesh", params2)]

    # DoF handler
    dof_handler_params = {"meshes": meshes}
    if model_type == ModelType.OnePhase:
      ic_params1 = {"mesh_name": meshes[0].name, "A": "0.2", "rho": "1", "u": "1", "p": "1"}
      ic_params2 = {"mesh_name": meshes[1].name, "A": "0.2", "rho": "1", "u": "1", "p": "1"}
      ics = list()
      ics.append(factory.createObject("InitialConditions1Phase", ic_params1))
      ics.append(factory.createObject("InitialConditions1Phase", ic_params2))
      dof_handler_class = "DoFHandler1Phase"
    else:
      ic_params1 = {"mesh_name": meshes[0].name, "A": "0.2", "vf1": "0.3",
        "rho1": "1", "u1": "1", "p1": "1", "rho2": "1", "u2": "1", "p2": "1"}
      ic_params2 = {"mesh_name": meshes[1].name, "A": "0.2", "vf1": "0.3",
        "rho1": "1", "u1": "1", "p1": "1", "rho2": "1", "u2": "1", "p2": "1"}
      ics = list()
      ics.append(factory.createObject("InitialConditions2Phase", ic_params1))
      ics.append(factory.createObject("InitialConditions2Phase", ic_params2))
      if model_type == ModelType.TwoPhaseNonInteracting:
        dof_handler_class = "DoFHandler2PhaseNonInteracting"
      elif model_type == ModelType.TwoPhase:
        dof_handler_class = "DoFHandler2Phase"
    dof_handler_params["ics"] = ics
    dof_handler = factory.createObject(dof_handler_class, dof_handler_params)

    # get indices of arhouA DoFs
    if model_type == ModelType.OnePhase:
      i_arhouA = [dof_handler.i(k, 2) for k in range(dof_handler.n_node)]
    elif model_type == ModelType.TwoPhaseNonInteracting:
      i_arhouA = [dof_handler.i(k, 2) for k in range(dof_handler.n_node)] \
        + [dof_handler.i(k, 5) for k in range(dof_handler.n_node)]
    else:
      i_arhouA = [dof_handler.i(k, 3) for k in range(dof_handler.n_node)] \
        + [dof_handler.i(k, 6) for k in range(dof_handler.n_node)]

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

    # update DoF handler with junction constraints
    dof_handler.updateWithJunctionConstraints([junction])
    n_dof = dof_handler.n_dof

    # compute base solution
    U = np.zeros(n_dof)
    U_old = np.zeros(n_dof)
    for i in range(n_dof):
      if use_zero_velocity and i in i_arhouA:
        U[i] = 0
        U_old[i] = 0
      else:
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
    for j in range(n_dof):
      # perturb solution
      U_perturbed = deepcopy(U)
      U_perturbed[j] += fd_eps

      # compute finite difference Jacobian
      r_perturbed = np.zeros(n_dof)
      J_perturbed = np.zeros(shape=(n_dof, n_dof))
      f(U_perturbed, U_old, r_perturbed, J_perturbed)
      for i in range(n_dof):
        J_fd[i,j] = (r_perturbed[i] - r[i]) / fd_eps

    # compute difference matrices
    abs_diffs = abs(J_hand_coded - J_fd)
    rel_diffs = computeRelativeDifferenceMatrix(J_hand_coded, J_fd)

    # print results
    print("\nJacobian, " + test_option + " contributions:")
    print("Hand-coded:")
    printMatrix(J_hand_coded)
    print("Finite-difference:")
    printMatrix(J_fd)
    print("Relative difference:")
    printRelativeMatrixDifference(rel_diffs, abs_diffs, 1e-1, 1e-3)

    matched = np.zeros((n_dof, n_dof), dtype=bool)
    for i in range(n_dof):
      for j in range(n_dof):
        matched[i,j] = abs_diffs[i,j] < self.abs_tol or rel_diffs[i,j] < self.rel_tol

    return matched
