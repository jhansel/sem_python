import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, VariableName

sys.path.append(base_dir + "src/closures")
from thermodynamic_functions import computeVolumeFraction, computeVelocity, \
  computeDensity, computeSpecificVolume, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy

sys.path.append(base_dir + "src/fem")
from FEValues import FEValues
from Quadrature import Quadrature

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class ExecutionerParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)

class Executioner(object):
  def __init__(self, params, model_type, ics, bcs, eos, interface_closures, gravity, dof_handler, mesh, nonlinear_solver_params, factory):
    self.model_type = model_type
    self.bcs = bcs
    self.eos = eos
    self.interface_closures = interface_closures
    self.gravity = gravity
    self.dof_handler = dof_handler
    self.mesh = mesh
    self.nonlinear_solver_params = nonlinear_solver_params
    self.quadrature = Quadrature()
    self.fe_values = FEValues(self.quadrature, dof_handler, mesh)
    self.factory = factory

    # initialize the solution
    self.U = np.zeros(self.dof_handler.n_dof)
    self.initializePhaseSolution(ics, 0)
    if self.model_type != ModelType.OnePhase:
      self.initializePhaseSolution(ics, 1)
    if self.model_type == ModelType.TwoPhase:
      self.initializeVolumeFractionSolution(ics)

    # create aux quantities
    self.aux1 = self.createIndependentPhaseAuxQuantities(0)
    if self.model_type != ModelType.OnePhase:
      self.aux2 = self.createIndependentPhaseAuxQuantities(1)
    if self.model_type == ModelType.TwoPhase:
      self.aux_2phase = self.createPhaseInteractionAuxQuantities()

    # create kernels
    self.kernels1 = self.createIndependentPhaseKernels(0)
    if self.model_type != ModelType.OnePhase:
      self.kernels2 = self.createIndependentPhaseKernels(1)

  def initializePhaseSolution(self, ics, phase):
    # get appropriate volume fraction function
    if self.model_type == ModelType.OnePhase:
      def initial_vf(x):
        return 1
    else:
      initial_vf1 = ics.vf1
      if phase == 0:
        initial_vf = initial_vf1
      else:
        def initial_vf(x):
          return 1 - initial_vf1(x)

    # get relevant IC functions
    initial_p = ics.p[phase]
    initial_u = ics.u[phase]
    if ics.specified_rho:
      initial_rho = ics.rho[phase]
    else:
      initial_T = ics.T[phase]

    # compute IC
    eos_phase = self.eos[phase]
    arho_index = self.dof_handler.variable_index[VariableName.ARho][phase]
    arhou_index = self.dof_handler.variable_index[VariableName.ARhoU][phase]
    arhoE_index = self.dof_handler.variable_index[VariableName.ARhoE][phase]
    for k in xrange(self.dof_handler.n_node):
      vf = initial_vf(self.mesh.x[k])
      p = initial_p(self.mesh.x[k])
      u = initial_u(self.mesh.x[k])
      if ics.specified_rho:
        rho = initial_rho(self.mesh.x[k])
      else:
        T = initial_T(self.mesh.x[k])
        rho = eos_phase.rho(p, T)
      e = eos_phase.e(1.0 / rho, p)[0]
      E = e + 0.5 * u * u
      self.U[self.dof_handler.i(k, arho_index)] = vf * rho
      self.U[self.dof_handler.i(k, arhou_index)] = vf * rho * u
      self.U[self.dof_handler.i(k, arhoE_index)] = vf * rho * E

  def initializeVolumeFractionSolution(self, ics):
    vf1_index = self.dof_handler.variable_index[VariableName.VF1]
    for k in xrange(self.dof_handler.n_node):
      self.U[self.dof_handler.i(k, vf1_index)] = ics.vf1(self.mesh.x[k])

  def createIndependentPhaseAuxQuantities(self, phase):
    # create list of aux quantities to create
    aux_names_phase = list()
    if phase == 0:
      if self.model_type == ModelType.OnePhase:
        aux_names_phase.append("VolumeFraction1Phase")
      else:
        aux_names_phase.append("VolumeFractionPhase1")
    else:
      aux_names_phase.append("VolumeFractionPhase2")
    aux_names_phase += ["Velocity", "SpecificTotalEnergy", "Density", \
                           "SpecificVolume", "SpecificInternalEnergy", "Pressure", "Temperature"]

    # create the aux quantities for this phase
    aux_list = list()
    for aux_name in aux_names_phase:
      params = {"phase": phase}
      if aux_name == "Pressure":
        params["p_function"] = self.eos[phase].p
      elif aux_name == "Temperature"
        params["T_function"] = self.eos[phase].T
      aux_list.append(self.factory.createObject(aux_name, params))

    return aux_list

  def createPhaseInteractionAuxQuantities(self):
    interaction_aux = list()

    aux_list = self.aux1 + self.aux2 + interaction_aux

    return aux_list

  def createIndependentPhaseKernels(self, phase):
    kernels = list()
    params = dict()
    params["phase"] = phase
    args = tuple([self.dof_handler])
    kernel_name_list = ["MassAdvection", "MomentumAdvection", "MomentumGravity", "EnergyAdvection", "EnergyGravity"]
    for kernel_name in kernel_name_list:
      kernels.append(self.factory.createObject(kernel_name, params, args))
    return kernels

  # computes the steady-state residual and Jacobian without applying strong BC
  def assembleSteadyStateSystemWithoutStrongBC(self, U):
    r = np.zeros(self.dof_handler.n_dof)
    J = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

    # volumetric terms
    self.addSteadyStateSystemPhase(U, self.kernels1, self.aux1, 0, r, J)
    if (self.model_type != ModelType.OnePhase):
      self.addSteadyStateSystemPhase(U, self.kernels2, self.aux2, 1, r, J)
    if (self.model_type == ModelType.TwoPhase):
      self.addSteadyStateSystemVolumeFraction(U, r, J)

    # weak boundary terms
    for bc in self.bcs:
      bc.applyWeakBC(U, r, J)

    return (r, J)

  # computes the full steady-state residual and Jacobian (strong BC applied)
  def assembleSteadyStateSystem(self, U):
    r, J = self.assembleSteadyStateSystemWithoutStrongBC(U)
    self.applyStrongBC(U, r, J)
    return (r, J)

  # applies strong BC
  def applyStrongBC(self, U, r, J):
    for bc in self.bcs:
      bc.applyStrongBC(U, r, J)

  # computes the steady-state residual and Jacobian for a phase
  def addSteadyStateSystemPhase(self, U, kernel_list, aux_list, phase, r, J):
    data = dict()
    der = dict()

    data["phi"] = self.fe_values.get_phi()
    data["g"] = self.gravity

    arho_name = "arho" + str(phase + 1)
    arhou_name = "arhou" + str(phase + 1)
    arhoE_name = "arhoE" + str(phase + 1)

    for elem in xrange(self.dof_handler.n_cell):
      r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
      J_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      data["grad_phi"] = self.fe_values.get_grad_phi(elem)
      data["JxW"] = self.fe_values.get_JxW(elem)

      # compute solution
      data["vf1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
      data[arho_name] = self.fe_values.computeLocalSolution(U, VariableName.ARho, phase, elem)
      data[arhou_name] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, phase, elem)
      data[arhoE_name] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, phase, elem)

      # compute auxiliary quantities
      for aux in aux_list:
        aux.compute(data, der)

      # compute the local residual and Jacobian
      for kernel in kernel_list:
        kernel.apply(data, der, r_cell, J_cell)

      # aggregate cell residual and matrix into global residual and matrix
      self.dof_handler.aggregateLocalVector(r, r_cell, elem)
      self.dof_handler.aggregateLocalMatrix(J, J_cell, elem)

  # computes the steady-state residual for the volume fraction equation
  def addSteadyStateSystemVolumeFraction(self, U, r, J):
    data = dict()
    der = dict()

    data["phi"] = self.fe_values.get_phi()
    data["g"] = self.gravity

    vf1_index = self.dof_handler.variable_index[VariableName.VF1]
    arho1_index = self.dof_handler.variable_index[VariableName.ARho][0]
    arhou1_index = self.dof_handler.variable_index[VariableName.ARhoU][0]
    arhoE1_index = self.dof_handler.variable_index[VariableName.ARhoE][0]
    arho2_index = self.dof_handler.variable_index[VariableName.ARho][1]
    arhou2_index = self.dof_handler.variable_index[VariableName.ARhoU][1]
    arhoE2_index = self.dof_handler.variable_index[VariableName.ARhoE][1]

    for elem in xrange(self.dof_handler.n_cell):
      r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
      J_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      data["grad_phi"] = self.fe_values.get_grad_phi(elem)
      data["JxW"] = self.fe_values.get_JxW(elem)

      # compute solution
      data["vf1"] = self.fe_values.computeLocalSolution(U, VariableName.VF1, 0, elem)
      data["arho1"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 0, elem)
      data["arhou1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 0, elem)
      data["arhoE1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 0, elem)
      data["arho2"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 1, elem)
      data["arhou2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 1, elem)
      data["arhoE2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 1, elem)

      # compute solution gradient
      data["dvf1_dx"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.VF1, 0, elem)

      # compute auxiliary quantities
      for aux in aux_list:
        aux.compute(data, der)

      # compute the local residual and Jacobian
      for kernel in kernel_list:
        kernel.apply(data, der, r_cell, J_cell)

      # compute the residual and Jacobian
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        i_vf1 = self.dof_handler.i(k_local, vf1_index)
        i_arhou1 = self.dof_handler.i(k_local, arhou1_index)
        i_arhoE1 = self.dof_handler.i(k_local, arhoE1_index)
        i_arhou2 = self.dof_handler.i(k_local, arhou2_index)
        i_arhoE2 = self.dof_handler.i(k_local, arhoE2_index)

        # compute local residual
        for q in xrange(self.quadrature.n_q):
          r_cell[i_vf1] += (uI[q] * dvf1_dx[q] - theta[q] * (p1[q] - p2[q])) * phi[k_local,q] * JxW[q]
          r_cell[i_arhou1] += - pI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]
          r_cell[i_arhoE1] += - pI[q] * uI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]
          r_cell[i_arhou2] += pI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]
          r_cell[i_arhoE2] += pI[q] * uI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]

        # compute local Jacobian
        for l_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
          j_vf1 = self.dof_handler.i(l_local, vf1_index)
          j_arho1 = self.dof_handler.i(l_local, arho1_index)
          j_arhou1 = self.dof_handler.i(l_local, arhou1_index)
          j_arhoE1 = self.dof_handler.i(l_local, arhoE1_index)
          j_arho2 = self.dof_handler.i(l_local, arho2_index)
          j_arhou2 = self.dof_handler.i(l_local, arhou2_index)
          j_arhoE2 = self.dof_handler.i(l_local, arhoE2_index)
          for q in xrange(self.quadrature.n_q):
            # volume fraction
            J_cell[i_vf1,j_vf1] += (uI[q] * grad_phi[l_local,q] \
              - dtheta_dvf1[q] * (p1[q] - p2[q]) * phi[l_local,q] \
              - theta[q] * (dp1_dvf1[q] - dp2_dvf1[q]) * phi[l_local,q]) * phi[k_local,q] * JxW[q]
            J_cell[i_vf1,j_arho1] += (duI_darho1[q] * dvf1_dx[q] - dtheta_darho1[q] * (p1[q] - p2[q]) \
              - theta[q] * dp1_darho1[q]) * phi[l_local,q] * phi[k_local,q] * JxW[q]
            J_cell[i_vf1,j_arhou1] += (duI_darhou1[q] * dvf1_dx[q] - dtheta_darhou1[q] * (p1[q] - p2[q]) \
              - theta[q] * dp1_darhou1[q]) * phi[l_local,q] * phi[k_local,q] * JxW[q]
            J_cell[i_vf1,j_arhoE1] += (- dtheta_darhoE1[q] * (p1[q] - p2[q]) \
              - theta[q] * dp1_darhoE1[q]) * phi[l_local,q] * phi[k_local,q] * JxW[q]
            J_cell[i_vf1,j_arho2] += (duI_darho2[q] * dvf1_dx[q] - dtheta_darho2[q] * (p1[q] - p2[q]) \
              + theta[q] * dp2_darho2[q]) * phi[l_local,q] * phi[k_local,q] * JxW[q]
            J_cell[i_vf1,j_arhou2] += (duI_darhou2[q] * dvf1_dx[q] - dtheta_darhou2[q] * (p1[q] - p2[q]) \
              + theta[q] * dp2_darhou2[q]) * phi[l_local,q] * phi[k_local,q] * JxW[q]
            J_cell[i_vf1,j_arhoE2] += (- dtheta_darhoE2[q] * (p1[q] - p2[q]) \
              + theta[q] * dp2_darhoE2[q]) * phi[l_local,q] * phi[k_local,q] * JxW[q]

            # momentum
            J[i_arhou1,j_vf1] += (- dpI_dvf1[q] * dvf1_dx[q] * phi[l_local,q] - pI[q] * grad_phi[l_local,q]) * phi[k_local,q] * JxW[q]
            J[i_arhou2,j_vf1] += (  dpI_dvf1[q] * dvf1_dx[q] * phi[l_local,q] + pI[q] * grad_phi[l_local,q]) * phi[k_local,q] * JxW[q]
            aux = dvf1_dx[q] * phi[l_local,q] * phi[k_local,q] * JxW[q]
            J[i_arhou1,j_arho1]  += - dpI_darho1[q]  * aux
            J[i_arhou2,j_arho1]  +=   dpI_darho1[q]  * aux
            J[i_arhou1,j_arhou1] += - dpI_darhou1[q] * aux
            J[i_arhou2,j_arhou1] +=   dpI_darhou1[q] * aux
            J[i_arhou1,j_arhoE1] += - dpI_darhoE1[q] * aux
            J[i_arhou2,j_arhoE1] +=   dpI_darhoE1[q] * aux
            J[i_arhou1,j_arho2]  += - dpI_darho2[q]  * aux
            J[i_arhou2,j_arho2]  +=   dpI_darho2[q]  * aux
            J[i_arhou1,j_arhou2] += - dpI_darhou2[q] * aux
            J[i_arhou2,j_arhou2] +=   dpI_darhou2[q] * aux
            J[i_arhou1,j_arhoE2] += - dpI_darhoE2[q] * aux
            J[i_arhou2,j_arhoE2] +=   dpI_darhoE2[q] * aux

            # energy
            J[i_arhoE1,j_vf1] += (- dpI_dvf1[q] * uI[q] * dvf1_dx[q] * phi[l_local,q] - pI[q] * uI[q] * grad_phi[l_local,q]) * phi[k_local,q] * JxW[q]
            J[i_arhoE2,j_vf1] += (  dpI_dvf1[q] * uI[q] * dvf1_dx[q] * phi[l_local,q] + pI[q] * uI[q] * grad_phi[l_local,q]) * phi[k_local,q] * JxW[q]
            J[i_arhoE1,j_arho1]  += - (dpI_darho1[q]  * uI[q] + pI[q] * duI_darho1[q])  * aux
            J[i_arhoE2,j_arho1]  +=   (dpI_darho1[q]  * uI[q] + pI[q] * duI_darho1[q])  * aux
            J[i_arhoE1,j_arhou1] += - (dpI_darhou1[q] * uI[q] + pI[q] * duI_darhou1[q]) * aux
            J[i_arhoE2,j_arhou1] +=   (dpI_darhou1[q] * uI[q] + pI[q] * duI_darhou1[q]) * aux
            J[i_arhoE1,j_arhoE1] += -  dpI_darhoE1[q] * uI[q]                           * aux
            J[i_arhoE2,j_arhoE1] +=    dpI_darhoE1[q] * uI[q]                           * aux
            J[i_arhoE1,j_arho2]  += - (dpI_darho2[q]  * uI[q] + pI[q] * duI_darho2[q])  * aux
            J[i_arhoE2,j_arho2]  +=   (dpI_darho2[q]  * uI[q] + pI[q] * duI_darho2[q])  * aux
            J[i_arhoE1,j_arhou2] += - (dpI_darhou2[q] * uI[q] + pI[q] * duI_darhou2[q]) * aux
            J[i_arhoE2,j_arhou2] +=   (dpI_darhou2[q] * uI[q] + pI[q] * duI_darhou2[q]) * aux
            J[i_arhoE1,j_arhoE2] += -  dpI_darhoE2[q] * uI[q]                           * aux
            J[i_arhoE2,j_arhoE2] +=    dpI_darhoE2[q] * uI[q]                           * aux

      # aggregate cell residual and matrix into global residual and matrix
      self.dof_handler.aggregateLocalVector(r, r_cell, elem)
      self.dof_handler.aggregateLocalMatrix(J, J_cell, elem)
