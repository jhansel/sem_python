from copy import deepcopy
import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from TestAux import TestAux, TestAuxParameters

sys.path.append(base_dir + "src/base")
from enums import ModelType, VariableName
from Factory import Factory

sys.path.append(base_dir + "src/fem")
from DoFHandler1Phase import DoFHandler1Phase
from DoFHandler2PhaseNonInteracting import DoFHandler2PhaseNonInteracting
from DoFHandler2Phase import DoFHandler2Phase
from FEValues import FEValues
from Quadrature import Quadrature

sys.path.append(base_dir + "src/mesh")
from UniformMesh import UniformMesh, UniformMeshParameters

sys.path.append(base_dir + "src/utilities")
from numeric_utilities import computeRelativeDifference

class KernelDerivativesTester(object):
  def __init__(self, verbose=False):
    self.verbose = verbose

  def checkDerivatives(self, kernel_name, model_type, phase, aux_dependencies, fd_eps=1e-8):
    self.model_type = model_type
    self.phase = phase

    # mesh
    params = UniformMeshParameters()
    params.set("n_cell", 1)
    mesh = UniformMesh(params)

    # DoF handler
    if self.model_type == ModelType.OnePhase:
      dof_handler = DoFHandler1Phase(mesh)
    elif self.model_type == ModelType.TwoPhaseNonInteracting:
      def vf1_initial(x):
        return 0.3
      dof_handler = DoFHandler2PhaseNonInteracting(mesh, vf1_initial)
    elif self.model_type == ModelType.TwoPhase:
      dof_handler = DoFHandler2Phase(mesh)

    # quadrature
    quadrature = Quadrature()

    # FE values
    self.fe_values = FEValues(quadrature, dof_handler, mesh)

    # factory
    factory = Factory()

    # kernel
    params = {"phase": phase}
    args = tuple([dof_handler])
    kernel = factory.createObject(kernel_name, params, args)

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
        for d,dependency in enumerate(aux_dependencies[aux_name]):
          coefs.append(a + 2.0 + d * 0.5)
        params.set("coefs", coefs)
        params.set("b", 1.0)
        aux_list.append(TestAux(params))

    # data
    data = dict()
    der = dict()
    self.elem = 0
    i = 0
    j = 1
    q = 0
    data["phi"] = self.fe_values.get_phi()
    data["grad_phi"] = self.fe_values.get_grad_phi(self.elem)
    data["JxW"] = self.fe_values.get_JxW(self.elem)

    # compute base solution
    U = np.zeros(dof_handler.n_dof)
    for k in xrange(dof_handler.n_dof):
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
    if self.verbose:
      for var_index in kernel.var_indices:
        var = dof_handler.variable_names[var_index]
        print "\nDerivative variable:", var
        print "  Hand-coded        =", J_hand_coded[var_index]
        print "  Finite difference =", J_fd[var_index]
        print "  Rel. difference   =", rel_diffs[var_index]

    # take the absolute value of the relative differences
    for x in rel_diffs:
      rel_diffs[x] = abs(rel_diffs[x])
    return rel_diffs

  def computeSolutionDependentData(self, U, data):
    data["vf1"] = self.fe_values.computeLocalVolumeFractionSolution(U, self.elem)
    data["grad_vf1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, self.elem)
    data["arho1"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 0, self.elem)
    data["arhou1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 0, self.elem)
    data["arhoE1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 0, self.elem)
    if self.model_type != ModelType.OnePhase:
      data["arho2"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 1, self.elem)
      data["arhou2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 1, self.elem)
      data["arhoE2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 1, self.elem)
