from Junction1Phase import Junction1Phase, Junction1PhaseParameters
from enums import ModelType
from error_utilities import error
from thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy

class CompressibleJunctionParameters(Junction1PhaseParameters):
  def __init__(self):
    Junction1PhaseParameters.__init__(self)
    self.registerBoolParameter("use_momentum_flux_balance", "Flag to use a momentum flux balance for last equation", False)

## Junction that uses compressible flow assumption
class CompressibleJunction(Junction1Phase):
  def __init__(self, params):
    Junction1Phase.__init__(self, params)
    if self.n_meshes != 2:
      error("CompressibleJunction is only implemented for connecting 2 meshes.")

    # determine if momentum flux balance is to be used for last equation;
    # otherwise, stagnation pressure condition will be used
    self.use_momentum_flux_balance = params.get("use_momentum_flux_balance")

    # get node indices
    self.k_left = self.node_indices[0]
    self.k_right = self.node_indices[1]
    self.k_left_adjacent = self.adjacent_node_indices[0]
    self.k_right_adjacent = self.adjacent_node_indices[1]

  def setDoFIndices(self):
    if (self.model_type == ModelType.TwoPhase):
      self.i_aA1L = self.dof_handler.i(self.k_left, self.aA1_index)
      self.i_aA1R = self.dof_handler.i(self.k_right, self.aA1_index)
      self.i_aA1LL = self.dof_handler.i(self.k_left_adjacent, self.aA1_index)
      self.i_aA1RR = self.dof_handler.i(self.k_right_adjacent, self.aA1_index)

    self.i_arhoL = self.dof_handler.i(self.k_left, self.arhoA_index)
    self.i_arhouAL = self.dof_handler.i(self.k_left, self.arhouA_index)
    self.i_arhoEAL = self.dof_handler.i(self.k_left, self.arhoEA_index)

    self.i_arhoR = self.dof_handler.i(self.k_right, self.arhoA_index)
    self.i_arhouAR = self.dof_handler.i(self.k_right, self.arhouA_index)
    self.i_arhoEAR = self.dof_handler.i(self.k_right, self.arhoEA_index)

    self.i_arhoLL = self.dof_handler.i(self.k_left_adjacent, self.arhoA_index)
    self.i_arhouALL = self.dof_handler.i(self.k_left_adjacent, self.arhouA_index)
    self.i_arhoEALL = self.dof_handler.i(self.k_left_adjacent, self.arhoEA_index)

    self.i_arhoRR = self.dof_handler.i(self.k_right_adjacent, self.arhoA_index)
    self.i_arhouARR = self.dof_handler.i(self.k_right_adjacent, self.arhouA_index)
    self.i_arhoEARR = self.dof_handler.i(self.k_right_adjacent, self.arhoEA_index)

  def applyWeaklyToNonlinearSystem(self, U, U_old, r, J):
    pass

  def applyStronglyToNonlinearSystem(self, U_new, U_old, r, J):
    n_values = 6
    L = 0
    R = 1
    L_old = 2
    R_old = 3
    LL_old = 4
    RR_old = 5

    k = [0] * n_values
    k[L] = self.k_left
    k[R] = self.k_right
    k[L_old] = self.k_left
    k[R_old] = self.k_right
    k[LL_old] = self.k_left_adjacent
    k[RR_old] = self.k_right_adjacent

    i_arhoA = [0] * n_values
    i_arhoA[L] = self.i_arhoL
    i_arhoA[R] = self.i_arhoR
    i_arhoA[L_old] = self.i_arhoL
    i_arhoA[R_old] = self.i_arhoR
    i_arhoA[LL_old] = self.i_arhoLL
    i_arhoA[RR_old] = self.i_arhoRR

    i_arhouA = [0] * n_values
    i_arhouA[L] = self.i_arhouAL
    i_arhouA[R] = self.i_arhouAR
    i_arhouA[L_old] = self.i_arhouAL
    i_arhouA[R_old] = self.i_arhouAR
    i_arhouA[LL_old] = self.i_arhouALL
    i_arhouA[RR_old] = self.i_arhouARR

    i_arhoEA = [0] * n_values
    i_arhoEA[L] = self.i_arhoEAL
    i_arhoEA[R] = self.i_arhoEAR
    i_arhoEA[L_old] = self.i_arhoEAL
    i_arhoEA[R_old] = self.i_arhoEAR
    i_arhoEA[LL_old] = self.i_arhoEALL
    i_arhoEA[RR_old] = self.i_arhoEARR

    U = [0] * n_values
    U[L] = U_new
    U[R] = U_new
    U[L_old] = U_old
    U[R_old] = U_old
    U[LL_old] = U_old
    U[RR_old] = U_old

    # initialize lists
    rho = [0] * n_values
    drho_daA1 = [0] * n_values
    drho_darhoA = [0] * n_values

    u = [0] * n_values
    du_darhoA = [0] * n_values
    du_darhouA = [0] * n_values

    e = [0] * n_values
    de_darhoA = [0] * n_values
    de_darhouA = [0] * n_values
    de_darhoEA = [0] * n_values

    p = [0] * n_values
    dp_daA1 = [0] * n_values
    dp_darhoA = [0] * n_values
    dp_darhouA = [0] * n_values
    dp_darhoEA = [0] * n_values

    c = [0] * n_values
    dc_daA1 = [0] * n_values
    dc_darhoA = [0] * n_values
    dc_darhouA = [0] * n_values
    dc_darhoEA = [0] * n_values

    if not self.use_momentum_flux_balance:
      p0 = [0] * n_values
      dp0_daA1 = [0] * n_values
      dp0_darhoA = [0] * n_values
      dp0_darhouA = [0] * n_values
      dp0_darhoEA = [0] * n_values

    # loop over subscript/superscript combinations
    for i in xrange(n_values):
      A = self.dof_handler.A[k[i]]
      aA1 = self.dof_handler.aA1(U[i], k[i])
      vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)
      arhoA  = U[i][i_arhoA[i]]
      arhouA = U[i][i_arhouA[i]]
      arhoEA = U[i][i_arhoEA[i]]

      rho[i], drho_dvf, drho_darhoA[i], _ = computeDensity(vf, arhoA, A)
      drho_daA1[i] = drho_dvf * dvf_daA1

      u[i], du_darhoA[i], du_darhouA[i] = computeVelocity(arhoA, arhouA)

      v, dv_drho = computeSpecificVolume(rho[i])
      dv_daA1 = dv_drho * drho_daA1[i]
      dv_darhoA = dv_drho * drho_darhoA[i]

      E, dE_darhoA, dE_darhoEA = computeSpecificTotalEnergy(arhoA, arhoEA)

      e[i], de_du, de_dE = computeSpecificInternalEnergy(u[i], E)
      de_darhoA[i] = de_du * du_darhoA[i] + de_dE * dE_darhoA
      de_darhouA[i] = de_du * du_darhouA[i]
      de_darhoEA[i] = de_dE * dE_darhoEA

      p[i], dp_dv, dp_de = self.eos.p(v, e[i])
      dp_daA1[i] = dp_dv * dv_daA1
      dp_darhoA[i] = dp_dv * dv_darhoA + dp_de * de_darhoA[i]
      dp_darhouA[i] = dp_de * de_darhouA[i]
      dp_darhoEA[i] = dp_de * de_darhoEA[i]

      c[i], dc_dv, dc_de = self.eos.c(v, e[i])
      dc_daA1[i] = dc_dv * dv_daA1
      dc_darhoA[i] = dc_dv * dv_darhoA + dc_de * de_darhoA[i]
      dc_darhouA[i] = dc_de * de_darhouA[i]
      dc_darhoEA[i] = dc_de * de_darhoEA[i]

      if not self.use_momentum_flux_balance:
        s, ds_dv, ds_de = self.eos.s(v, e[i])
        ds_daA1 = ds_dv * dv_daA1
        ds_darhoA = ds_dv * dv_darhoA + ds_de * de_darhoA[i]
        ds_darhouA = ds_de * de_darhouA[i]
        ds_darhoEA = ds_de * de_darhoEA[i]

        h, dh_de, dh_dp, dh_drho = computeSpecificEnthalpy(e[i], p[i], rho[i])
        dh_daA1 = dh_dp * dp_daA1[i]
        dh_darhoA = dh_de * de_darhoA[i] + dh_dp * dp_darhoA[i] + dh_drho * drho_darhoA[i]
        dh_darhouA = dh_de * de_darhouA[i] + dh_dp * dp_darhouA[i]
        dh_darhoEA = dh_de * de_darhoEA[i] + dh_dp * dp_darhoEA[i]

        h0 = h + 0.5 * u[i]**2
        dh0_daA1 = dh_daA1
        dh0_darhoA = dh_darhoA + u[i] * du_darhoA[i]
        dh0_darhouA = dh_darhouA + u[i] * du_darhouA[i]
        dh0_darhoEA = dh_darhoEA

        p0[i], dp0_dh0, dp0_ds = self.eos.p_from_h_s(h0, s)
        dp0_daA1[i] = dp0_dh0 * dh0_daA1 + dp0_ds * ds_daA1
        dp0_darhoA[i] = dp0_dh0 * dh0_darhoA + dp0_ds * ds_darhoA
        dp0_darhouA[i] = dp0_dh0 * dh0_darhouA + dp0_ds * ds_darhouA
        dp0_darhoEA[i] = dp0_dh0 * dh0_darhoEA + dp0_ds * ds_darhoEA

    # compute old average quantities
    rhoL = 0.5 * (rho[L_old] + rho[LL_old])
    rhoR = 0.5 * (rho[R_old] + rho[RR_old])
    uL = 0.5 * (u[L_old] + u[LL_old])
    uR = 0.5 * (u[R_old] + u[RR_old])
    pL = 0.5 * (p[L_old] + p[LL_old])
    pR = 0.5 * (p[R_old] + p[RR_old])
    cL = 0.5 * (c[L_old] + c[LL_old])
    cR = 0.5 * (c[R_old] + c[RR_old])

    # residual indices
    i1 = self.i_arhoL
    i2 = self.i_arhouAL
    i3 = self.i_arhoEAL
    i4 = self.i_arhoR
    i5 = self.i_arhouAR
    i6 = self.i_arhoEAR

    # reset Jacobian rows
    J[i1,:] = 0
    J[i2,:] = 0
    J[i3,:] = 0
    J[i4,:] = 0
    J[i5,:] = 0
    J[i6,:] = 0

    # compute residuals and Jacobians
    r[i1] = p[L] - pL + rhoL * (cL - uL) * (u[L] - uL)
    J_i1_vf1L = dp_daA1[L]
    J[i1,self.i_arhoL] = dp_darhoA[L] + rhoL * (cL - uL) * du_darhoA[L]
    J[i1,self.i_arhouAL] = dp_darhouA[L] + rhoL * (cL - uL) * du_darhouA[L]
    J[i1,self.i_arhoEAL] = dp_darhoEA[L]

    r[i2] = p[R] - pR - rhoR * (cR - uR) * (u[R] - uR)
    J_i2_vf1R = dp_daA1[R]
    J[i2,self.i_arhoR] = dp_darhoA[R] - rhoR * (cR - uR) * du_darhoA[R]
    J[i2,self.i_arhouAR] = dp_darhouA[R] - rhoR * (cR - uR) * du_darhouA[R]
    J[i2,self.i_arhoEAR] = dp_darhoEA[R]

    if u[L] >= 0:
      if u[R] < 0:
        error("Assumption violated: Both velocity conditions were true.")
      r[i3] = rho[L] - rhoL - (p[L] - pL) / cL**2
      J_i3_vf1L = drho_daA1[L] - dp_daA1[L] / cL**2
      J_i3_vf1R = 0
      J[i3,self.i_arhoL] = drho_darhoA[L] - dp_darhoA[L] / cL**2
      J[i3,self.i_arhouAL] = - dp_darhouA[L] / cL**2
      J[i3,self.i_arhoEAL] = - dp_darhoEA[L] / cL**2
    elif u[R] < 0:
      r[i3] = rho[R] - rhoR - (p[R] - pR) / cR**2
      J_i3_vf1R = drho_daA1[R] - dp_daA1[R] / cR**2
      J_i3_vf1L = 0
      J[i3,self.i_arhoR] = drho_darhoA[R] - dp_darhoA[R] / cR**2
      J[i3,self.i_arhouAR] = - dp_darhouA[R] / cR**2
      J[i3,self.i_arhoEAR] = - dp_darhoEA[R] / cR**2
    else:
      error("Assumption violated: Neither velocity condition was true.")

    r[i4] = rho[L] * u[L] - rho[R] * u[R]
    J_i4_vf1L = drho_daA1[L] * u[L]
    J[i4,self.i_arhoL] = drho_darhoA[L] * u[L] + rho[L] * du_darhoA[L]
    J[i4,self.i_arhouAL] = rho[L] * du_darhouA[L]
    J_i4_vf1R = -drho_daA1[R] * u[R]
    J[i4,self.i_arhoR] = -(drho_darhoA[R] * u[R] + rho[R] * du_darhoA[R])
    J[i4,self.i_arhouAR] = -rho[R] * du_darhouA[R]

    r[i5] = e[L] + p[L] / rho[L] + 0.5 * u[L]**2 - (e[R] + p[R] / rho[R] + 0.5 * u[R]**2)
    J_i5_vf1L = dp_daA1[L] / rho[L] - p[L] / rho[L]**2 * drho_daA1[L]
    J[i5,self.i_arhoL] = de_darhoA[L] + dp_darhoA[L] / rho[L] - p[L] / rho[L]**2 * drho_darhoA[L] + u[L] * du_darhoA[L]
    J[i5,self.i_arhouAL] = de_darhouA[L] + dp_darhouA[L] / rho[L] + u[L] * du_darhouA[L]
    J[i5,self.i_arhoEAL] = de_darhoEA[L] + dp_darhoEA[L] / rho[L]
    J_i5_vf1R = -(dp_daA1[R] / rho[R] - p[R] / rho[R]**2 * drho_daA1[R])
    J[i5,self.i_arhoR] = -(de_darhoA[R] + dp_darhoA[R] / rho[R] - p[R] / rho[R]**2 * drho_darhoA[R] + u[R] * du_darhoA[R])
    J[i5,self.i_arhouAR] = -(de_darhouA[R] + dp_darhouA[R] / rho[R] + u[R] * du_darhouA[R])
    J[i5,self.i_arhoEAR] = -(de_darhoEA[R] + dp_darhoEA[R] / rho[R])

    if self.use_momentum_flux_balance:
      r[i6] = rho[L] * u[L]**2 + p[L] - (rho[R] * u[R]**2 + p[R])
      J_i6_vf1L = drho_daA1[L] * u[L]**2 + dp_daA1[L]
      J[i6,self.i_arhoL] = drho_darhoA[L] * u[L]**2 + rho[L] * 2.0 * u[L] * du_darhoA[L] + dp_darhoA[L]
      J[i6,self.i_arhouAL] = rho[L] * 2.0 * u[L] * du_darhouA[L] + dp_darhouA[L]
      J[i6,self.i_arhoEAL] = dp_darhoEA[L]
      J_i6_vf1R = -(drho_daA1[R] * u[R]**2 + dp_daA1[R])
      J[i6,self.i_arhoR] = -(drho_darhoA[R] * u[R]**2 + rho[R] * 2.0 * u[R] * du_darhoA[R] + dp_darhoA[R])
      J[i6,self.i_arhouAR] = -(rho[R] * 2.0 * u[R] * du_darhouA[R] + dp_darhouA[R])
      J[i6,self.i_arhoEAR] = -dp_darhoEA[R]
    else:
      r[i6] = p0[L] - p0[R]
      J_i6_vf1L = dp0_daA1[L]
      J[i6,self.i_arhoL] = dp0_darhoA[L]
      J[i6,self.i_arhouAL] = dp0_darhouA[L]
      J[i6,self.i_arhoEAL] = dp0_darhoEA[L]
      J_i6_vf1R = -dp0_daA1[R]
      J[i6,self.i_arhoR] = -dp0_darhoA[R]
      J[i6,self.i_arhouAR] = -dp0_darhouA[R]
      J[i6,self.i_arhoEAR] = -dp0_darhoEA[R]

    if self.model_type == ModelType.TwoPhase:
      J[i1,self.i_aA1L] = J_i1_vf1L
      J[i2,self.i_aA1R] = J_i2_vf1R
      J[i3,self.i_aA1L] = J_i3_vf1L
      J[i3,self.i_aA1R] = J_i3_vf1R
      J[i4,self.i_aA1L] = J_i4_vf1L
      J[i4,self.i_aA1R] = J_i4_vf1R
      J[i5,self.i_aA1L] = J_i5_vf1L
      J[i5,self.i_aA1R] = J_i5_vf1R
      J[i6,self.i_aA1L] = J_i6_vf1L
      J[i6,self.i_aA1R] = J_i6_vf1R

  def applyStronglyToLinearSystemMatrix(self, A):
    pass

  def applyStronglyToLinearSystemRHSVector(self, U_old, b):
    pass
