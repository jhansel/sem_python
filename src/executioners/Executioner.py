import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, VariableName

sys.path.append(base_dir + "src/closures")
from thermodynamic_functions import computeVelocity, computeDensity, \
  computeSpecificVolume, computeSpecificTotalEnergy, computeSpecificInternalEnergy

sys.path.append(base_dir + "src/fem")
from FEValues import FEValues
from Quadrature import Quadrature

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class ExecutionerParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)

class Executioner(object):
  def __init__(self, params, model_type, ics, bcs, eos, interface_closures, gravity, dof_handler, mesh, nonlinear_solver_params):
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

    # initialize the solution
    self.U = np.zeros(self.dof_handler.n_dof)
    if self.model_type == ModelType.OnePhase:
      self.initializeOnePhaseSolution(ics)
    else:
      self.initializeTwoPhaseSolution(ics, 0)
      self.initializeTwoPhaseSolution(ics, 1)
      if self.model_type == ModelType.TwoPhase:
        self.initializeVolumeFractionSolution(ics)

  def initializeOnePhaseSolution(self, ics):
    phase = 0
    p0 = ics.p0[phase]
    if ics.specified_rho:
      rho0 = ics.rho0[phase]
    else:
      T0 = ics.T0[phase]
    u0 = ics.u0[phase]
    eos_phase = self.eos[phase]
    for k in xrange(self.dof_handler.n_node):
      p = p0(self.mesh.x[k])
      u = u0(self.mesh.x[k])
      if ics.specified_rho:
        rho = rho0(self.mesh.x[k])
      else:
        T = T0(self.mesh.x[k])
        rho = eos_phase.rho(p, T)
      e = eos_phase.e(1.0 / rho, p)[0]
      E = e + 0.5 * u * u
      self.U[self.dof_handler.i(k, VariableName.ARho, phase)] = rho
      self.U[self.dof_handler.i(k, VariableName.ARhoU, phase)] = rho * u
      self.U[self.dof_handler.i(k, VariableName.ARhoE, phase)] = rho * E

  def initializeTwoPhaseSolution(self, ics, phase):
    vf0 = ics.vf0
    p0 = ics.p0[phase]
    if ics.specified_rho:
      rho0 = ics.rho0[phase]
    else:
      T0 = ics.T0[phase]
    u0 = ics.u0[phase]
    eos_phase = self.eos[phase]
    for k in xrange(self.dof_handler.n_node):
      vf1 = vf0(self.mesh.x[k])
      if phase == 0:
        vf = vf1
      else:
        vf = 1 - vf1
      p = p0(self.mesh.x[k])
      u = u0(self.mesh.x[k])
      if ics.specified_rho:
        rho = rho0(self.mesh.x[k])
      else:
        T = T0(self.mesh.x[k])
        rho = eos_phase.rho(p, T)
      e = eos_phase.e(1.0 / rho, p)[0]
      E = e + 0.5 * u * u
      self.U[self.dof_handler.i(k, VariableName.ARho, phase)] = vf * rho
      self.U[self.dof_handler.i(k, VariableName.ARhoU, phase)] = vf * rho * u
      self.U[self.dof_handler.i(k, VariableName.ARhoE, phase)] = vf * rho * E

  def initializeVolumeFractionSolution(self, ics):
    for k in xrange(self.dof_handler.n_node):
      self.U[self.dof_handler.i(k, VariableName.VF1)] = ics.vf0(self.mesh.x[k])

  # computes the steady-state residual and Jacobian without applying strong BC
  def assembleSteadyStateSystemWithoutStrongBC(self, U):
    r = np.zeros(self.dof_handler.n_dof)
    J = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

    # volumetric terms
    self.addSteadyStateSystemPhase(U, 0, r, J)
    if (self.model_type != ModelType.OnePhase):
      self.addSteadyStateSystemPhase(U, 1, r, J)
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
  def addSteadyStateSystemPhase(self, U, phase, r, J):
    phi = self.fe_values.get_phi()

    for elem in xrange(self.dof_handler.n_cell):
      r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
      J_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      grad_phi = self.fe_values.get_grad_phi(elem)
      JxW = self.fe_values.get_JxW(elem)

      # compute solution
      vf, dvf_dvf1 = self.fe_values.computeLocalVolumeFractionSolution(U, phase, elem)
      arho = self.fe_values.computeLocalSolution(U, VariableName.ARho, phase, elem)
      arhou = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, phase, elem)
      arhoE = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, phase, elem)

      # compute auxiliary quantities
      u, du_darho, du_darhou = computeVelocity(arho, arhou)

      rho, drho_dvf, drho_darho = computeDensity(vf, arho)
      drho_dvf1 = drho_dvf * dvf_dvf1

      v, dv_drho = computeSpecificVolume(rho)
      dv_dvf1 = dv_drho * drho_dvf1
      dv_darho = dv_drho * drho_darho

      E, dE_darho, dE_darhoE = computeSpecificTotalEnergy(arho, arhoE)

      e, de_du, de_dE = computeSpecificInternalEnergy(u, E)
      de_darho = de_du * du_darho + de_dE * dE_darho
      de_darhou = de_du * du_darhou
      de_darhoE = de_dE * dE_darhoE

      p, dp_dv, dp_de = self.eos[phase].p(v, e)
      dp_dvf1 = dp_dv * dv_dvf1
      dp_darho = dp_dv * dv_darho + dp_de * de_darho
      dp_darhou = dp_de * de_darhou
      dp_darhoE = dp_de * de_darhoE

      # compute the residual and Jacobian
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        i_arho = self.dof_handler.i(k_local, VariableName.ARho, phase)
        i_arhou = self.dof_handler.i(k_local, VariableName.ARhoU, phase)
        i_arhoE = self.dof_handler.i(k_local, VariableName.ARhoE, phase)

        # compute local residual
        for q in xrange(self.quadrature.n_q):
          r_cell[i_arho] += - arhou[q] * grad_phi[k_local,q] * JxW[q]
          r_cell[i_arhou] += (- (arhou[q] * u[q] + vf[q] * p[q]) * grad_phi[k_local,q] - arho[q] * self.gravity * phi[k_local,q]) * JxW[q]
          r_cell[i_arhoE] += (- u[q] * (arhoE[q] + vf[q] * p[q]) * grad_phi[k_local,q] - arhou[q] * self.gravity * phi[k_local,q]) * JxW[q]

        # compute local Jacobian
        for l_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
          if (self.model_type == ModelType.TwoPhase):
            j_vf1 = self.dof_handler.i(l_local, VariableName.VF1)
          j_arho = self.dof_handler.i(l_local, VariableName.ARho, phase)
          j_arhou = self.dof_handler.i(l_local, VariableName.ARhoU, phase)
          j_arhoE = self.dof_handler.i(l_local, VariableName.ARhoE, phase)
          for q in xrange(self.quadrature.n_q):
            # mass
            J_cell[i_arho,j_arhou] += - phi[l_local,q] * grad_phi[k_local,q] * JxW[q]

            # momentum
            if (self.model_type == ModelType.TwoPhase):
              J_cell[i_arhou,j_vf1] += - (dvf_dvf1[q] * p[q] + vf[q] * dp_dvf1[q]) * phi[l_local,q] * grad_phi[k_local,q] * JxW[q]
            J_cell[i_arhou,j_arho] += (- (arhou[q] * du_darho[q] + vf[q] * dp_darho[q]) * phi[l_local,q] * grad_phi[k_local,q] \
              - self.gravity * phi[l_local,q] * phi[k_local,q]) * JxW[q]
            J_cell[i_arhou,j_arhou] += - (u[q] + arhou[q] * du_darhou[q] + vf[q] * dp_darhou[q]) * phi[l_local,q] * grad_phi[k_local,q] * JxW[q]
            J_cell[i_arhou,j_arhoE] += - vf[q] * dp_darhoE[q] * phi[l_local,q] * grad_phi[k_local,q] * JxW[q]

            # energy
            if (self.model_type == ModelType.TwoPhase):
              J_cell[i_arhoE,j_vf1] += - u[q] * (vf[q] * dp_dvf1[q] + dvf_dvf1[q] * p[q]) * phi[l_local,q] * grad_phi[k_local,q] * JxW[q]
            J_cell[i_arhoE,j_arho] += - (u[q] * (vf[q] * dp_darho[q]) + du_darho[q] * (arhoE[q] + vf[q] * p[q])) * phi[l_local,q] * grad_phi[k_local,q] * JxW[q]
            J_cell[i_arhoE,j_arhou] += (- (u[q] * (vf[q] * dp_darhou[q]) + du_darhou[q] * (arhoE[q] + vf[q] * p[q])) * phi[l_local,q] * grad_phi[k_local,q] \
              - self.gravity * phi[l_local,q] * phi[k_local,q]) * JxW[q]
            J_cell[i_arhoE,j_arhoE] += - u[q] * (1 + vf[q] * dp_darhoE[q]) * phi[l_local,q] * grad_phi[k_local,q] * JxW[q]

      # aggregate cell residual and matrix into global residual and matrix
      self.dof_handler.aggregateLocalVector(r, r_cell, elem)
      self.dof_handler.aggregateLocalMatrix(J, J_cell, elem)

  # computes the steady-state residual for the volume fraction equation
  def addSteadyStateSystemVolumeFraction(self, U, r, J):
    phi = self.fe_values.get_phi()

    for elem in xrange(self.dof_handler.n_cell):
      r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
      J_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      grad_phi = self.fe_values.get_grad_phi(elem)
      JxW = self.fe_values.get_JxW(elem)

      phase1 = 0
      phase2 = 1
      eos1 = self.eos[phase1]
      eos2 = self.eos[phase2]

      # compute solution
      vf1 = self.fe_values.computeLocalSolution(U, VariableName.VF1, phase1, elem)
      arho1 = self.fe_values.computeLocalSolution(U, VariableName.ARho, phase1, elem)
      arhou1 = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, phase1, elem)
      arhoE1 = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, phase1, elem)
      arho2 = self.fe_values.computeLocalSolution(U, VariableName.ARho, phase2, elem)
      arhou2 = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, phase2, elem)
      arhoE2 = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, phase2, elem)

      # compute solution gradient
      dvf1_dx = self.fe_values.computeLocalSolutionGradient(U, VariableName.VF1, phase1, elem)

      # compute auxiliary quantities
      vf2 = 1 - vf1
      dvf2_dvf1 = 0 * vf2 - 1

      rho1, drho1_dvf1, drho1_darho1 = computeDensity(vf1, arho1)
      rho2, drho2_dvf2, drho2_darho2 = computeDensity(vf2, arho2)
      drho2_dvf1 = drho2_dvf2 * dvf2_dvf1

      v1, dv1_drho1 = computeSpecificVolume(rho1)
      v2, dv2_drho2 = computeSpecificVolume(rho2)
      dv1_dvf1 = dv1_drho1 * drho1_dvf1
      dv2_dvf1 = dv2_drho2 * drho2_dvf1
      dv1_darho1 = dv1_drho1 * drho1_darho1
      dv2_darho2 = dv2_drho2 * drho2_darho2

      u1, du1_darho1, du1_darhou1 = computeVelocity(arho1, arhou1)
      u2, du2_darho2, du2_darhou2 = computeVelocity(arho2, arhou2)

      E1, dE1_darho1, dE1_darhoE1 = computeSpecificTotalEnergy(arho1, arhoE1)
      E2, dE2_darho2, dE2_darhoE2 = computeSpecificTotalEnergy(arho2, arhoE2)

      e1, de1_du1, de1_dE1 = computeSpecificInternalEnergy(u1, E1)
      e2, de2_du2, de2_dE2 = computeSpecificInternalEnergy(u2, E2)
      de1_darho1 = de1_du1 * du1_darho1 + de1_dE1 * dE1_darho1
      de2_darho2 = de2_du2 * du2_darho2 + de2_dE2 * dE2_darho2
      de1_darhou1 = de1_du1 * du1_darhou1
      de2_darhou2 = de2_du2 * du2_darhou2
      de1_darhoE1 = de1_dE1 * dE1_darhoE1
      de2_darhoE2 = de2_dE2 * dE2_darhoE2

      p1, dp1_dv1, dp1_de1 = eos1.p(v1, e1)
      p2, dp2_dv2, dp2_de2 = eos2.p(v2, e2)
      dp1_dvf1 = dp1_dv1 * dv1_dvf1
      dp2_dvf1 = dp2_dv2 * dv2_dvf1
      dp1_darho1 = dp1_dv1 * dv1_darho1 + dp1_de1 * de1_darho1
      dp2_darho2 = dp2_dv2 * dv2_darho2 + dp2_de2 * de2_darho2
      dp1_darhou1 = dp1_de1 * de1_darhou1
      dp2_darhou2 = dp2_de2 * de2_darhou2
      dp1_darhoE1 = dp1_de1 * de1_darhoE1
      dp2_darhoE2 = dp2_de2 * de2_darhoE2

      T1, dT1_dv1, dT1_de1 = eos1.T(v1, e1)
      T2, dT2_dv2, dT2_de2 = eos2.T(v2, e2)
      dT1_dvf1 = dT1_dv1 * dv1_dvf1
      dT2_dvf1 = dT2_dv2 * dv2_dvf1
      dT1_darho1 = dT1_dv1 * dv1_darho1 + dT1_de1 * de1_darho1
      dT2_darho2 = dT2_dv2 * dv2_darho2 + dT2_de2 * de2_darho2
      dT1_darhou1 = dT1_de1 * de1_darhou1
      dT2_darhou2 = dT2_de2 * de2_darhou2
      dT1_darhoE1 = dT1_de1 * de1_darhoE1
      dT2_darhoE2 = dT2_de2 * de2_darhoE2

      beta, dbeta_darho1, dbeta_darho2 = self.interface_closures.computeBeta(arho1, arho2)

      mu, dmu_dT1, dmu_dT2, dmu_dbeta = self.interface_closures.computeMu(T1, T2, beta)
      dmu_dvf1 = dmu_dT1 * dT1_dvf1 + dmu_dT2 * dT2_dvf1
      dmu_darho1 = dmu_dT1 * dT1_darho1 + dmu_dbeta * dbeta_darho1
      dmu_darho2 = dmu_dT2 * dT2_darho2 + dmu_dbeta * dbeta_darho2
      dmu_darhou1 = dmu_dT1 * dT1_darhou1
      dmu_darhou2 = dmu_dT2 * dT2_darhou2
      dmu_darhoE1 = dmu_dT1 * dT1_darhoE1
      dmu_darhoE2 = dmu_dT2 * dT2_darhoE2

      uI, duI_du1, duI_du2, duI_dbeta = self.interface_closures.computeInterfaceVelocity(u1, u2, beta)
      duI_darho1 = duI_du1 * du1_darho1 + duI_dbeta * dbeta_darho1
      duI_darho2 = duI_du2 * du2_darho2 + duI_dbeta * dbeta_darho2
      duI_darhou1 = duI_du1 * du1_darhou1
      duI_darhou2 = duI_du2 * du2_darhou2

      pI, dpI_dp1, dpI_dp2, dpI_dmu = self.interface_closures.computeInterfacePressure(p1, p2, mu)
      dpI_dvf1 = dpI_dp1 * dp1_dvf1 + dpI_dp2 * dp2_dvf1 + dpI_dmu * dmu_dvf1
      dpI_darho1 = dpI_dp1 * dp1_darho1 + dpI_dmu * dmu_darho1
      dpI_darho2 = dpI_dp2 * dp2_darho2 + dpI_dmu * dmu_darho2
      dpI_darhou1 = dpI_dp1 * dp1_darhou1 + dpI_dmu * dmu_darhou1
      dpI_darhou2 = dpI_dp2 * dp2_darhou2 + dpI_dmu * dmu_darhou2
      dpI_darhoE1 = dpI_dp1 * dp1_darhoE1 + dpI_dmu * dmu_darhoE1
      dpI_darhoE2 = dpI_dp2 * dp2_darhoE2 + dpI_dmu * dmu_darhoE2

      theta, pdtheta_pdvf1, dtheta_dp1, dtheta_dp2 = self.interface_closures.computeTheta(vf1, p1, p2)
      dtheta_dvf1 = pdtheta_pdvf1 + dtheta_dp1 * dp1_dvf1 + dtheta_dp2 * dp2_dvf1
      dtheta_darho1 = dtheta_dp1 * dp1_darho1
      dtheta_darho2 = dtheta_dp2 * dp2_darho2
      dtheta_darhou1 = dtheta_dp1 * dp1_darhou1
      dtheta_darhou2 = dtheta_dp2 * dp2_darhou2
      dtheta_darhoE1 = dtheta_dp1 * dp1_darhoE1
      dtheta_darhoE2 = dtheta_dp2 * dp2_darhoE2

      # compute the residual and Jacobian
      for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
        i_vf1 = self.dof_handler.i(k_local, VariableName.VF1, phase1)
        i_arhou1 = self.dof_handler.i(k_local, VariableName.ARhoU, phase1)
        i_arhoE1 = self.dof_handler.i(k_local, VariableName.ARhoE, phase1)
        i_arhou2 = self.dof_handler.i(k_local, VariableName.ARhoU, phase2)
        i_arhoE2 = self.dof_handler.i(k_local, VariableName.ARhoE, phase2)

        # compute local residual
        for q in xrange(self.quadrature.n_q):
          r_cell[i_vf1] += (uI[q] * dvf1_dx[q] - theta[q] * (p1[q] - p2[q])) * phi[k_local,q] * JxW[q]
          r_cell[i_arhou1] += - pI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]
          r_cell[i_arhoE1] += - pI[q] * uI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]
          r_cell[i_arhou2] += pI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]
          r_cell[i_arhoE2] += pI[q] * uI[q] * dvf1_dx[q] * phi[k_local,q] * JxW[q]

        # compute local Jacobian
        for l_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
          j_vf1 = self.dof_handler.i(l_local, VariableName.VF1, phase1)
          j_arho1 = self.dof_handler.i(l_local, VariableName.ARho, phase1)
          j_arhou1 = self.dof_handler.i(l_local, VariableName.ARhoU, phase1)
          j_arhoE1 = self.dof_handler.i(l_local, VariableName.ARhoE, phase1)
          j_arho2 = self.dof_handler.i(l_local, VariableName.ARho, phase2)
          j_arhou2 = self.dof_handler.i(l_local, VariableName.ARhoU, phase2)
          j_arhoE2 = self.dof_handler.i(l_local, VariableName.ARhoE, phase2)
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
