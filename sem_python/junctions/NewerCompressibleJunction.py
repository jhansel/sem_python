from .Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..base.enums import ModelType
from ..utilities.error_utilities import error
from ..closures.thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy

class NewerCompressibleJunctionParameters(Junction1PhaseParameters):
  def __init__(self):
    Junction1PhaseParameters.__init__(self)
    self.registerFloatListParameter("loss_coefficients", "List of loss coefficients for each mesh")

## Junction that uses compressible flow assumption
class NewerCompressibleJunction(Junction1Phase):
  def __init__(self, params):
    Junction1Phase.__init__(self, params)

    # get loss coefficients, if any
    if params.has("loss_coefficients"):
      self.loss_coefficients = params.get("loss_coefficients")
      if len(self.loss_coefficients) != self.n_meshes:
        error("The list parameters 'loss_coefficients' and 'mesh_names' must have the same size.")
    else:
      self.loss_coefficients = [0] * self.n_meshes

    self.n_constraints += 1

    self.initializeVariableVectors()

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
    self.computeJunctionStagnationEnthalpy()

    # add the boundary fluxes; characteristic theory is used to decide how
    # many pieces of information are supplied to inlets vs. outlets
    for i in xrange(self.n_meshes):
      if self.u[i] * self.nx[i] > 0: # inlets to the junction; outlet BC
        self.addJunctionInletFlux(i, U_new, r, J)
      else: # outlets from the junction; inlet BC
        self.addJunctionOutletFlux(i, U_new, r, J)

    # the normal boundary fluxes are computed only because the mass fluxes are
    # needed by the mass flux balance constraint
    self.computeFluxes(U_new)

  ## Computes junction stagnation enthalpy and its derivatives
  def computeJunctionStagnationEnthalpy(self):
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
      if self.u[i] * self.nx[i] > 0: # mesh serves as inlet
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

  def addJunctionOutletFlux(self, i, U, r, J):
    # get the junction entropy
    sJ = U[self.i_sJ]

    # extract solution data
    k = self.node_indices[i]

    A = self.dof_handler.A[k]
    arhouA = U[self.i_arhouA[i]]

    nx = self.nx[i]

    vf = self.vf[i]
    u = self.u[i]

    H = self.h0J
    dH_darhoA  = self.dh0J_darhoA[i]
    dH_darhouA = self.dh0J_darhouA[i]
    dH_darhoEA = self.dh0J_darhoEA[i]

    p0J, dp0J_dh0J, dp0J_dsJ = self.eos.p_from_h_s(self.h0J, sJ)
    dp0J_darhoA  = dp0J_dh0J * self.dh0J_darhoA[i]
    dp0J_darhouA = dp0J_dh0J * self.dh0J_darhouA[i]
    dp0J_darhoEA = dp0J_dh0J * self.dh0J_darhoEA[i]

    p0 = p0J - self.loss_coefficients[i] * (p0J - self.p[i])
    dp0_dp0J = 1 - self.loss_coefficients[i]
    dp0_dp = self.loss_coefficients[i]
    dp0_daA1    = dp0_dp * self.dp_daA1[i]
    dp0_darhoA  = dp0_dp0J * dp0J_darhoA  + dp0_dp * self.dp_darhoA[i]
    dp0_darhouA = dp0_dp0J * dp0J_darhouA + dp0_dp * self.dp_darhouA[i]
    dp0_darhoEA = dp0_dp0J * dp0J_darhoEA + dp0_dp * self.dp_darhoEA[i]
    dp0_dsJ = dp0_dp0J * dp0J_dsJ

    s, ds_dH, ds_dp0 = self.eos.s_from_h_p(H, p0)
    ds_daA1 = ds_dp0 * dp0_daA1
    ds_darhoA  = ds_dH * dH_darhoA  + ds_dp0 * dp0_darhoA
    ds_darhouA = ds_dH * dH_darhouA + ds_dp0 * dp0_darhouA
    ds_darhoEA = ds_dH * dH_darhoEA + ds_dp0 * dp0_darhoEA
    ds_dsJ = ds_dp0 * dp0_dsJ

    h = H - 0.5 * u**2
    dh_darhoA  = dH_darhoA  - u * self.du_darhoA[i]
    dh_darhouA = dH_darhouA - u * self.du_darhouA[i]
    dh_darhoEA = dH_darhoEA

    p, dp_dh, dp_ds = self.eos.p_from_h_s(h, s)
    dp_daA1 = dp_ds * ds_daA1
    dp_darhoA  = dp_dh * dh_darhoA  + dp_ds * ds_darhoA
    dp_darhouA = dp_dh * dh_darhouA + dp_ds * ds_darhouA
    dp_darhoEA = dp_dh * dh_darhoEA + dp_ds * ds_darhoEA
    dp_dsJ = dp_ds * ds_dsJ

    rho, drho_dp, drho_ds = self.eos.rho_from_p_s(p, s)
    drho_daA1    = drho_dp * dp_daA1    + drho_ds * ds_daA1
    drho_darhoA  = drho_dp * dp_darhoA  + drho_ds * ds_darhoA
    drho_darhouA = drho_dp * dp_darhouA + drho_ds * ds_darhouA
    drho_darhoEA = drho_dp * dp_darhoEA + drho_ds * ds_darhoEA
    drho_dsJ     = drho_dp * dp_dsJ     + drho_ds * ds_dsJ

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

    r[i_energy] += vf * u * rho * H * A * nx
    J[i_energy][self.i_sJ] += vf * u * drho_dsJ * H * A * nx
    J[i_energy][i_mass] += vf * (self.du_darhoA[i] * rho * H \
      + u * (drho_darhoA * H + rho * dH_darhoA)) * A * nx
    J[i_energy][i_momentum] += vf * (self.du_darhouA[i] * rho * H \
      + u * (drho_darhouA * H + rho * dH_darhouA)) * A * nx
    J[i_energy][i_energy] += vf * u * (drho_darhoEA * H + rho * dH_darhoEA) * A * nx

    # add contributions from other inlets through h0J
    for j in xrange(self.n_meshes):
      if j != i:
        dH_darhoAj  = self.dh0J_darhoA[j]
        dH_darhouAj = self.dh0J_darhouA[j]
        dH_darhoEAj = self.dh0J_darhoEA[j]

        dp0J_darhoAj  = dp0J_dh0J * dH_darhoAj
        dp0J_darhouAj = dp0J_dh0J * dH_darhouAj
        dp0J_darhoEAj = dp0J_dh0J * dH_darhoEAj

        dp0_darhoAj  = dp0_dp0J * dp0J_darhoAj
        dp0_darhouAj = dp0_dp0J * dp0J_darhouAj
        dp0_darhoEAj = dp0_dp0J * dp0J_darhoEAj

        ds_darhoAj  = ds_dH * dH_darhoAj  + ds_dp0 * dp0_darhoAj
        ds_darhouAj = ds_dH * dH_darhouAj + ds_dp0 * dp0_darhouAj
        ds_darhoEAj = ds_dH * dH_darhoEAj + ds_dp0 * dp0_darhoEAj

        dh_darhoAj  = dH_darhoAj
        dh_darhouAj = dH_darhouAj
        dh_darhoEAj = dH_darhoEAj

        dp_darhoAj  = dp_dh * dh_darhoAj  + dp_ds * ds_darhoAj
        dp_darhouAj = dp_dh * dh_darhouAj + dp_ds * ds_darhouAj
        dp_darhoEAj = dp_dh * dh_darhoEAj + dp_ds * ds_darhoEAj

        drho_darhoAj  = drho_dp * dp_darhoAj  + drho_ds * ds_darhoAj
        drho_darhouAj = drho_dp * dp_darhouAj + drho_ds * ds_darhouAj
        drho_darhoEAj = drho_dp * dp_darhoEAj + drho_ds * ds_darhoEAj

        j_mass = self.i_arhoA[j]
        j_momentum = self.i_arhouA[j]
        j_energy = self.i_arhoEA[j]

        J[i_mass][j_mass]     += vf * drho_darhoAj  * u * A * nx
        J[i_mass][j_momentum] += vf * drho_darhouAj * u * A * nx
        J[i_mass][j_energy]   += vf * drho_darhoEAj * u * A * nx

        J[i_momentum][j_mass]     += vf * (drho_darhoAj  * u**2 + dp_darhoAj)  * A * nx
        J[i_momentum][j_momentum] += vf * (drho_darhouAj * u**2 + dp_darhouAj) * A * nx
        J[i_momentum][j_energy]   += vf * (drho_darhoEAj * u**2 + dp_darhoEAj) * A * nx

        J[i_energy][j_mass]     += vf * u * (drho_darhoAj  * H + rho * dH_darhoAj)  * A * nx
        J[i_energy][j_momentum] += vf * u * (drho_darhouAj * H + rho * dH_darhouAj) * A * nx
        J[i_energy][j_energy]   += vf * u * (drho_darhoEAj * H + rho * dH_darhoEAj) * A * nx

  def addJunctionInletFlux(self, i, U, r, J):
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
