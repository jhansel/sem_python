from copy import deepcopy
import numpy as np

from sem_python.aux.TestAux import TestAux, TestAuxParameters
from sem_python.base.enums import ModelType, VariableName
from sem_python.base.Factory import Factory
from sem_python.utilities.numeric_utilities import computeRelativeDifference

class KernelDerivativesTester(object):
  def checkDerivatives(self, kernel_name, model_type, phase, aux_dependencies, \
                       aux_gradients=list(), kernel_params=dict(), fd_eps=1e-8):
    self.model_type = model_type
    self.phase = phase

    # factory
    factory = Factory()

    # mesh
    params = {"n_cell": 1}
    meshes = [factory.createObject("UniformMesh", params)]

    # DoF handler
    dof_handler_params = {"meshes": meshes}
    if self.model_type == ModelType.OnePhase:
      ic_params = {"mesh_name": meshes[0].name, "A": "2", "rho": "1", "u": "1", "p": "1"}
      ics = [factory.createObject("InitialConditions1Phase", ic_params)]
      dof_handler_class = "DoFHandler1Phase"
    else:
      ic_params = {"mesh_name": meshes[0].name, "A": "2", "vf1": "0.3", "rho1": "1",
                   "u1": "1", "p1": "1", "rho2": "1", "u2": "1", "p2": "1"}
      ics = [factory.createObject("InitialConditions2Phase", ic_params)]
      if self.model_type == ModelType.TwoPhaseNonInteracting:
        dof_handler_class = "DoFHandler2PhaseNonInteracting"
      elif self.model_type == ModelType.TwoPhase:
        dof_handler_class = "DoFHandler2Phase"
    dof_handler_params["ics"] = ics
    dof_handler = factory.createObject(dof_handler_class, dof_handler_params)

    # quadrature
    quadrature_params = {}
    quadrature = factory.createObject("Quadrature", quadrature_params)

    # FE values
    fe_values_params = {"quadrature": quadrature, "dof_handler": dof_handler, "meshes": meshes}
    self.fe_values = factory.createObject("FEValues", fe_values_params)

    # kernel
    kernel_params["phase"] = phase
    kernel_params["dof_handler"] = dof_handler
    kernel = factory.createObject(kernel_name, kernel_params)

    # aux
    aux_list = list()
    for a,aux_name in enumerate(aux_dependencies):
      # "vf1" is a special case of aux because its name is also the name of
      # a solution variable (in 2-phase); therefore one needs to make sure that
      # it uses its own "identity" aux instead of the generic test aux
      if aux_name == "vf1":
        params = {"phase": 0}
        aux_list.append(factory.createObject("VolumeFractionPhase1", params))
      else:
        params = TestAuxParameters()
        params.set("var", aux_name)
        params.set("other_vars", aux_dependencies[aux_name])
        coefs = list()
        for d, _ in enumerate(aux_dependencies[aux_name]):
          coefs.append(a + 2.0 + d * 0.5)
        params.set("coefs", coefs)
        params.set("b", 1.0)
        aux_list.append(TestAux(params))

    # add the aux derivatives
    for aux_name in aux_gradients:
      params = {"aux": aux_name, "variable_names": aux_dependencies[aux_name]}
      aux_list.append(factory.createObject("AuxGradient", params))

    # data
    data = dict()
    aux_names = [aux.name for aux in aux_list]
    der = dof_handler.initializeDerivativeData(aux_names)
    for aux_name in aux_dependencies:
      for dep in aux_dependencies[aux_name]:
        der[aux_name][dep] = 0
    self.elem = 0
    i = 0
    j = 1
    q = 0
    data["htc_wall"] = 0.2
    data["T_wall"] = 1.3
    data["P_heat"] = 0.6
    data["grad_A"] = 0.3
    data["phi"] = self.fe_values.get_phi()
    data["grad_phi"] = self.fe_values.get_grad_phi(self.elem)
    data["JxW"] = self.fe_values.get_JxW(self.elem)

    # compute base solution
    U = np.zeros(dof_handler.n_dof)
    for k in range(dof_handler.n_dof):
      U[k] = k + 1.0
    self.computeSolutionDependentData(U, data)
    for aux in aux_list:
      aux.compute(data, der)

    # base calculation
    r = kernel.computeResidual(data, i)[q]
    J_hand_coded = dict()
    for var_index in kernel.var_indices:
      J_hand_coded[var_index] = kernel.computeJacobian(data, der, var_index, i, j)[q]

    # finite difference Jacobians
    rel_diffs = dict()
    J_fd = dict()
    for var_index in kernel.var_indices:
      # perturb solution and recompute aux
      U_perturbed = deepcopy(U)
      j_global = dof_handler.i(j, var_index)
      U_perturbed[j_global] += fd_eps
      self.computeSolutionDependentData(U_perturbed, data)
      for aux in aux_list:
        aux.compute(data, der)

      # compute finite difference Jacobian
      r_perturbed = kernel.computeResidual(data, i)[q]
      J_fd[var_index] = (r_perturbed - r) / fd_eps
      rel_diffs[var_index] = computeRelativeDifference(J_hand_coded[var_index], J_fd[var_index])

    # print results
    for var_index in kernel.var_indices:
      var = dof_handler.variable_names[var_index]
      print("\nDerivative variable:", var)
      print("  Hand-coded        =", J_hand_coded[var_index])
      print("  Finite difference =", J_fd[var_index])
      print("  Rel. difference   =", rel_diffs[var_index])

    # take the absolute value of the relative differences
    for x in rel_diffs:
      rel_diffs[x] = abs(rel_diffs[x])
    return rel_diffs

  def computeSolutionDependentData(self, U, data):
    data["A"] = self.fe_values.computeLocalArea(self.elem)
    data["aA1"] = self.fe_values.computeLocalVolumeFractionSolution(U, self.elem)
    data["arhoA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 0, self.elem)
    data["arhouA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 0, self.elem)
    data["arhoEA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 0, self.elem)
    data["grad_aA1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, self.elem)
    data["grad_arhoA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoA, 0, self.elem)
    data["grad_arhouA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoUA, 0, self.elem)
    data["grad_arhoEA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoEA, 0, self.elem)
    if self.model_type != ModelType.OnePhase:
      data["arhoA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 1, self.elem)
      data["arhouA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 1, self.elem)
      data["arhoEA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 1, self.elem)
      data["grad_arhoA2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoA, 1, self.elem)
      data["grad_arhouA2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoUA, 1, self.elem)
      data["grad_arhoEA2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoEA, 1, self.elem)
