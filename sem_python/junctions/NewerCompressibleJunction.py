from Junction1Phase import Junction1Phase, Junction1PhaseParameters
from enums import ModelType
from error_utilities import error
from thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy

class NewerCompressibleJunctionParameters(Junction1PhaseParameters):
  def __init__(self):
    Junction1PhaseParameters.__init__(self)

## Junction that uses compressible flow assumption
class NewerCompressibleJunction(Junction1Phase):
  def __init__(self, params):
    Junction1Phase.__init__(self, params)

    self.n_constraints += 1

  def setDoFIndices(self):
    Junction1Phase.setDoFIndices(self)

    self.i_constraint_mass = self.i_constraint[0]
    self.i_sJ = self.i_constraint_mass

  def initializeConstraintVariables(self, U):
    # initialize the junction h0 and s with the average of all IC values for these quantities
    s_sum = 0
    for i in xrange(self.n_meshes):
      k = self.node_indices[i]

      aA1 = self.dof_handler.aA1(U, k)
      arhoA = U[self.i_arhoA[i]]
      arhouA = U[self.i_arhouA[i]]
      arhoEA = U[self.i_arhoEA[i]]

      vf, _ = computeVolumeFraction(aA1, self.A[i], self.phase, self.model_type)
      rho, _, _, _ = computeDensity(vf, arhoA, self.A[i])
      v, _ = computeSpecificVolume(rho)
      u, _, _ = computeVelocity(arhoA, arhouA)
      E, _, _ = computeSpecificTotalEnergy(arhoA, arhoEA)
      e, _, _ = computeSpecificInternalEnergy(u, E)
      s, _, _ = self.eos.s(v, e)

      s_sum += s

    U[self.i_sJ] = s_sum / self.n_meshes

  def applyWeaklyToNonlinearSystem(self, U_new, U_old, r, J):
    self.computeFlowQuantities(U_new)
    self.computeJunctionStagnationEnthalpy(U_old)

    # add the boundary fluxes; characteristic theory is used to decide how
    # many pieces of information are supplied to inlets vs. outlets
    for i in xrange(self.n_meshes):
      arhouA_old = U_old[self.i_arhouA[i]]
      if arhouA_old * self.nx[i] > 0: # inlets to the junction; outlet BC
        self.addOutletBC(i, U_new, r, J)
      else: # outlets from the junction; inlet BC
        self.addInletBC(i, U_new, r, J)

    # the normal boundary fluxes are computed only because the mass fluxes are
    # needed by the mass flux balance constraint
    self.computeFluxes(U_new)

  ## Computes junction stagnation enthalpy and its derivatives
  def computeJunctionStagnationEnthalpy(self, U_old):
    # compute junction stagnation enthalpy h0J and its derivatives
    num = 0
    den = 0
    num_alt = 0
    den_alt = 0
    dnum_darhoA  = [0] * self.n_meshes
    dnum_darhouA = [0] * self.n_meshes
    dnum_darhoEA = [0] * self.n_meshes
    dden_darhoA  = [0] * self.n_meshes
    dden_darhouA = [0] * self.n_meshes
    dden_darhoEA = [0] * self.n_meshes
    dnum_alt_darhoA  = [0] * self.n_meshes
    dnum_alt_darhouA = [0] * self.n_meshes
    dnum_alt_darhoEA = [0] * self.n_meshes
    for i in xrange(self.n_meshes):
      arhouA_old = U_old[self.i_arhouA[i]]
      if arhouA_old * self.nx[i] > 0: # mesh serves as inlet
        num += self.rho[i] * self.u[i] * self.h0[i] * self.A[i] * self.nx[i]
        den += self.rho[i] * self.u[i] * self.A[i] * self.nx[i]
        dnum_darhoA[i] = (self.drho_darhoA[i] * self.u[i] * self.h0[i] \
          + self.rho[i] * self.du_darhoA[i] * self.h0[i] \
          + self.rho[i] * self.u[i] * self.dh0_darhoA[i]) * self.A[i] * self.nx[i]
        dnum_darhouA[i] = self.rho[i] * (self.du_darhouA[i] * self.h0[i] \
          + self.u[i] * self.dh0_darhouA[i]) * self.A[i] * self.nx[i]
        dnum_darhoEA[i] = self.rho[i] * self.u[i] * self.dh0_darhoEA[i] * self.A[i] * self.nx[i]
        dden_darhoA[i] = (self.drho_darhoA[i] * self.u[i] + self.rho[i] * self.du_darhoA[i]) * self.A[i] * self.nx[i]
        dden_darhouA[i] = self.rho[i] * self.du_darhouA[i] * self.A[i] * self.nx[i]
        dden_darhoA[i] = 0
      else: # add to separate numerator / denominator in case there is little or no flow into junction
        num_alt += self.h0[i]
        dnum_alt_darhoA[i]  = self.dh0_darhoA[i]
        dnum_alt_darhouA[i] = self.dh0_darhouA[i]
        dnum_alt_darhoEA[i] = self.dh0_darhoEA[i]
        den_alt += 1

    # finish computing h0J and its derivatives
    self.dh0J_darhoA  = [0] * self.n_meshes
    self.dh0J_darhouA = [0] * self.n_meshes
    self.dh0J_darhoEA = [0] * self.n_meshes
    if den > 1e-10:
      # use normal definition
      self.h0J = num / den
      for i in xrange(self.n_meshes):
        self.dh0J_darhoA[i]  = dnum_darhoA[i]  / den - num / den**2 * dden_darhoA[i]
        self.dh0J_darhouA[i] = dnum_darhouA[i] / den - num / den**2 * dden_darhouA[i]
        self.dh0J_darhoEA[i] = dnum_darhoEA[i] / den - num / den**2 * dden_darhoEA[i]
    else: # little or no flow into the junction
      # use alternate definition
      self.h0J = num_alt / den_alt
      for i in xrange(self.n_meshes):
        self.dh0J_darhoA[i]  = dnum_alt_darhoA[i]  / den_alt
        self.dh0J_darhouA[i] = dnum_alt_darhouA[i] / den_alt
        self.dh0J_darhoEA[i] = dnum_alt_darhoEA[i] / den_alt

  def addInletBC(self, i, U, r, J):
    # get the junction entropy
    sJ = U[self.i_sJ]

    # extract solution data
    k = self.node_indices[i]

    A = self.dof_handler.A[k]
    arhouA = U[self.i_arhouA[i]]

    nx = self.nx[i]

    vf = self.vf[i]
    u = self.u[i]

    h = self.h0J - 0.5 * u**2
    dh_darhoA  = self.dh0J_darhoA[i]  -u * self.du_darhoA[i]
    dh_darhouA = self.dh0J_darhouA[i] -u * self.du_darhouA[i]
    dh_darhoEA = self.dh0J_darhoEA[i]

    p, dp_dh, dp_dsJ = self.eos.p_from_h_s(h, sJ)
    dp_darhoA  = dp_dh * dh_darhoA
    dp_darhouA = dp_dh * dh_darhouA
    dp_darhoEA = dp_dh * dh_darhoEA

    rho, drho_dp, drho_dsJ_partial = self.eos.rho_from_p_s(p, sJ)
    drho_dsJ = drho_dsJ_partial + drho_dp * dp_dsJ
    drho_darhoA  = drho_dp * dp_darhoA
    drho_darhouA = drho_dp * dp_darhouA
    drho_darhoEA = drho_dp * dp_darhoEA

    v, dv_drho = computeSpecificVolume(rho)
    dv_dsJ = dv_drho * drho_dsJ
    dv_darhoA  = dv_drho * drho_darhoA
    dv_darhouA = dv_drho * drho_darhouA
    dv_darhoEA = dv_drho * drho_darhoEA

    e, de_dv, de_dp = self.eos.e(v, p)
    de_dsJ = de_dv * dv_dsJ + de_dp * dp_dsJ
    de_darhoA  = de_dv * dv_darhoA  + de_dp * dp_darhoA
    de_darhouA = de_dv * dv_darhouA + de_dp * dp_darhouA
    de_darhoEA = de_dv * dv_darhoEA + de_dp * dp_darhoEA

    E = e + 0.5 * u**2
    dE_dsJ = de_dsJ
    dE_darhoA  = de_darhoA  + u * self.du_darhoA[i]
    dE_darhouA = de_darhouA + u * self.du_darhouA[i]
    dE_darhoEA = de_darhoEA

    # compute the fluxes
    i_mass = self.i_arhoA[i]
    i_momentum = self.i_arhouA[i]
    i_energy = self.i_arhoEA[i]

    r[i_mass] += vf * rho * u * A * nx
    J[i_mass][self.i_sJ] += vf * drho_dsJ * u * A * nx
    J[i_mass][i_mass] += vf * (drho_darhoA * u + rho * self.du_darhoA[i]) * A * nx
    J[i_mass][i_momentum] += vf * (drho_darhouA * u + rho * self.du_darhouA[i]) * A * nx
    J[i_mass][i_energy] += vf * drho_darhoEA * u * A * nx

    r[i_momentum] += vf * (rho * u**2 + p) * A * nx
    J[i_momentum][self.i_sJ] += vf * (drho_dsJ * u**2 + dp_dsJ) * A * nx
    J[i_momentum][i_mass] += vf * (drho_darhoA * u**2 + rho * 2 * u * self.du_darhoA[i] + dp_darhoA) * A * nx
    J[i_momentum][i_momentum] += vf * (drho_darhouA * u**2 + rho * 2 * u * self.du_darhouA[i] + dp_darhouA) * A * nx
    J[i_momentum][i_energy] += vf * (drho_darhoEA * u**2 + dp_darhoEA) * A * nx

    r[i_energy] += vf * u * (rho * E + p) * A * nx
    J[i_energy][self.i_sJ] += vf * u * (drho_dsJ * E + rho * dE_dsJ + dp_dsJ) * A * nx
    J[i_energy][i_mass] += vf * (self.du_darhoA[i] * (rho * E + p) \
      + u * (drho_darhoA * E + rho * dE_darhoA + dp_darhoA)) * A * nx
    J[i_energy][i_momentum] += vf * (self.du_darhouA[i] * (rho * E + p) \
      + u * (drho_darhouA * E + rho * dE_darhouA + dp_darhouA)) * A * nx
    J[i_energy][i_energy] += vf * u * (drho_darhoEA * E + rho * dE_darhoEA + dp_darhoEA) * A * nx

    # add contributions from other inlets through h0J
    for j in xrange(self.n_meshes):
      if j != i:
        dh_darhoAj  = self.dh0J_darhoA[j]
        dh_darhouAj = self.dh0J_darhouA[j]
        dh_darhoEAj = self.dh0J_darhoEA[j]

        dp_darhoAj  = dp_dh * dh_darhoAj
        dp_darhouAj = dp_dh * dh_darhouAj
        dp_darhoEAj = dp_dh * dh_darhoEAj

        drho_darhoAj  = drho_dp * dp_darhoAj
        drho_darhouAj = drho_dp * dp_darhouAj
        drho_darhoEAj = drho_dp * dp_darhoEAj

        dv_darhoAj  = dv_drho * drho_darhoAj
        dv_darhouAj = dv_drho * drho_darhouAj
        dv_darhoEAj = dv_drho * drho_darhoEAj

        de_darhoAj  = de_dv * dv_darhoAj  + de_dp * dp_darhoAj
        de_darhouAj = de_dv * dv_darhouAj + de_dp * dp_darhouAj
        de_darhoEAj = de_dv * dv_darhoEAj + de_dp * dp_darhoEAj

        dE_darhoAj  = de_darhoAj
        dE_darhouAj = de_darhouAj
        dE_darhoEAj = de_darhoEAj

        j_mass = self.i_arhoA[j]
        j_momentum = self.i_arhouA[j]
        j_energy = self.i_arhoEA[j]

        J[i_mass][j_mass]     += vf * drho_darhoAj  * u * A * nx
        J[i_mass][j_momentum] += vf * drho_darhouAj * u * A * nx
        J[i_mass][j_energy]   += vf * drho_darhoEAj * u * A * nx

        J[i_momentum][j_mass]     += vf * (drho_darhoAj  * u**2 + dp_darhoAj)  * A * nx
        J[i_momentum][j_momentum] += vf * (drho_darhouAj * u**2 + dp_darhouAj) * A * nx
        J[i_momentum][j_energy]   += vf * (drho_darhoEAj * u**2 + dp_darhoEAj) * A * nx

        J[i_energy][j_mass]     += vf * u * (drho_darhoAj  * E + rho * dE_darhoAj  + dp_darhoAj)  * A * nx
        J[i_energy][j_momentum] += vf * u * (drho_darhouAj * E + rho * dE_darhouAj + dp_darhouAj) * A * nx
        J[i_energy][j_energy]   += vf * u * (drho_darhoEAj * E + rho * dE_darhoEAj + dp_darhoEAj) * A * nx

  def addOutletBC(self, i, U, r, J):
    # get the junction entropy
    sJ = U[self.i_sJ]

    # extract solution data
    k = self.node_indices[i]

    A = self.dof_handler.A[k]
    arhouA = U[self.i_arhouA[i]]

    nx = self.nx[i]

    vf = self.vf[i]
    rho = self.rho[i]
    v = self.v[i]
    u = self.u[i]

    h = self.h0J - 0.5 * u**2
    dh_darhoA  = self.dh0J_darhoA[i]  -u * self.du_darhoA[i]
    dh_darhouA = self.dh0J_darhouA[i] -u * self.du_darhouA[i]
    dh_darhoEA = self.dh0J_darhoEA[i]

    p, dp_dh, dp_dsJ = self.eos.p_from_h_s(h, sJ)
    dp_darhoA  = dp_dh * dh_darhoA
    dp_darhouA = dp_dh * dh_darhouA
    dp_darhoEA = dp_dh * dh_darhoEA

    e, de_dv, de_dp = self.eos.e(v, p)
    de_dsJ = de_dp * dp_dsJ
    de_darhoA  = de_dv * self.dv_darhoA[i] + de_dp * dp_darhoA
    de_darhouA = de_dp * dp_darhouA
    de_darhoEA = de_dp * dp_darhoEA

    E = e + 0.5 * u**2
    dE_dsJ = de_dsJ
    dE_darhoA  = de_darhoA  + u * self.du_darhoA[i]
    dE_darhouA = de_darhouA + u * self.du_darhouA[i]
    dE_darhoEA = de_darhoEA

    # compute the fluxes
    i_mass = self.i_arhoA[i]
    i_momentum = self.i_arhouA[i]
    i_energy = self.i_arhoEA[i]

    r[i_mass] += arhouA * nx
    J[i_mass][i_momentum] += nx

    r[i_momentum] += (arhouA * u + vf * p * A) * nx
    J[i_momentum][self.i_sJ] += vf * dp_dsJ * A * nx
    J[i_momentum][i_mass] += (arhouA * self.du_darhoA[i] + vf * dp_darhoA * A) * nx
    J[i_momentum][i_momentum] += (arhouA * self.du_darhouA[i] + u + vf * dp_darhouA * A) * nx
    J[i_momentum][i_energy] += vf * dp_darhoEA * A * nx

    r[i_energy] += vf * u * (rho * E + p) * A * nx
    J[i_energy][self.i_sJ] += vf * u * (rho * dE_dsJ + dp_dsJ) * A * nx
    J[i_energy][i_mass] += vf * (self.du_darhoA[i] * (rho * E + p) \
      + u * (self.drho_darhoA[i] * E + rho * dE_darhoA + dp_darhoA)) * A * nx
    J[i_energy][i_momentum] += vf * (self.du_darhouA[i] * (rho * E + p) \
      + u * (rho * dE_darhouA + dp_darhouA)) * A * nx
    J[i_energy][i_energy] += vf * u * (rho * dE_darhoEA + dp_darhoEA) * A * nx

    # add contributions from other inlets through h0J
    for j in xrange(self.n_meshes):
      if j != i:
        dh_darhoAj  = self.dh0J_darhoA[j]
        dh_darhouAj = self.dh0J_darhouA[j]
        dh_darhoEAj = self.dh0J_darhoEA[j]

        dp_darhoAj  = dp_dh * dh_darhoAj
        dp_darhouAj = dp_dh * dh_darhouAj
        dp_darhoEAj = dp_dh * dh_darhoEAj

        de_darhoAj  = de_dp * dp_darhoAj
        de_darhouAj = de_dp * dp_darhouAj
        de_darhoEAj = de_dp * dp_darhoEAj

        dE_darhoAj  = de_darhoAj
        dE_darhouAj = de_darhouAj
        dE_darhoEAj = de_darhoEAj

        j_mass = self.i_arhoA[j]
        j_momentum = self.i_arhouA[j]
        j_energy = self.i_arhoEA[j]

        J[i_momentum][j_mass]     += vf * dp_darhoAj  * A * nx
        J[i_momentum][j_momentum] += vf * dp_darhouAj * A * nx
        J[i_momentum][j_energy]   += vf * dp_darhoEAj * A * nx

        J[i_energy][j_mass]     += vf * u * (rho * dE_darhoAj  + dp_darhoAj)  * A * nx
        J[i_energy][j_momentum] += vf * u * (rho * dE_darhouAj + dp_darhouAj) * A * nx
        J[i_energy][j_energy]   += vf * u * (rho * dE_darhoEAj + dp_darhoEAj) * A * nx

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    r[self.i_constraint_mass] = sum(self.f_mass)
    J[self.i_constraint_mass,:] = 0
    for n in xrange(self.n_meshes):
      J[self.i_constraint_mass,self.i_arhouA[n]] = self.df_mass_darhouA[n]

  def applyStronglyToLinearSystemMatrix(self, A):
    pass

  def applyStronglyToLinearSystemRHSVector(self, U_old, b):
    pass
