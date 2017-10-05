from .Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..base.enums import ModelType
from ..utilities.error_utilities import error
from ..closures.thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy, addKineticEnergy, \
  subtractKineticEnergy

class NewestCompressibleJunctionParameters(Junction1PhaseParameters):
  def __init__(self):
    Junction1PhaseParameters.__init__(self)

## Junction that uses compressible flow assumption
class NewestCompressibleJunction(Junction1Phase):
  def __init__(self, params):
    Junction1Phase.__init__(self, params)

    self.n_constraints += 3

    self.f_mass = [0] * self.n_meshes
    self.df_mass_daA1 = [0] * self.n_meshes
    self.df_mass_darhoA = [0] * self.n_meshes
    self.df_mass_darhouA = [0] * self.n_meshes
    self.df_mass_dHJ = [0] * self.n_meshes
    self.df_mass_dsJ = [0] * self.n_meshes

    self.f_energy = [0] * self.n_meshes
    self.df_energy_daA1 = [0] * self.n_meshes
    self.df_energy_darhoA = [0] * self.n_meshes
    self.df_energy_darhouA = [0] * self.n_meshes
    self.df_energy_dpJ = [0] * self.n_meshes
    self.df_energy_dHJ = [0] * self.n_meshes
    self.df_energy_dsJ = [0] * self.n_meshes

    self.f_entropy = [0] * self.n_meshes
    self.df_entropy_daA1 = [0] * self.n_meshes
    self.df_entropy_darhoA = [0] * self.n_meshes
    self.df_entropy_darhouA = [0] * self.n_meshes
    self.df_entropy_dpJ = [0] * self.n_meshes
    self.df_entropy_dHJ = [0] * self.n_meshes
    self.df_entropy_dsJ = [0] * self.n_meshes

    self.dH_avg_darhoA  = [0] * self.n_meshes
    self.dH_avg_darhouA = [0] * self.n_meshes
    self.dH_avg_darhoEA = [0] * self.n_meshes

    self.ds_avg_darhoA  = [0] * self.n_meshes
    self.ds_avg_darhouA = [0] * self.n_meshes
    self.ds_avg_darhoEA = [0] * self.n_meshes

  def setDoFIndices(self):
    Junction1Phase.setDoFIndices(self)

    self.i_constraint_mass = self.i_constraint[0]
    self.i_constraint_energy = self.i_constraint[1]
    self.i_constraint_entropy = self.i_constraint[2]
    self.i_pJ = self.i_constraint_mass
    self.i_HJ = self.i_constraint_energy
    self.i_sJ = self.i_constraint_entropy

  def initializeConstraintVariables(self, U):
    # initialize the junction variables with the average of all IC values for these quantities
    p_sum = 0
    H_sum = 0
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
      p, _, _ = self.eos.p(v, e)
      h, _, _, _ = computeSpecificEnthalpy(e, p, rho)
      H, _, _ = addKineticEnergy(h, u)
      s, _, _ = self.eos.s(v, e)

      p_sum += p
      H_sum += H
      s_sum += s

    U[self.i_pJ] = p_sum / self.n_meshes
    U[self.i_HJ] = H_sum / self.n_meshes
    U[self.i_sJ] = s_sum / self.n_meshes

  def addJunctionInletFlux(self, i, U, r, J):
    k = self.node_indices[i]

    A = self.A[i]
    nx = self.nx[i]
    aA1 = self.dof_handler.aA1(U, k)
    arhoA = U[self.i_arhoA[i]]
    arhouA = U[self.i_arhouA[i]]
    pJ = U[self.i_pJ]

    vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

    rho, drho_dvf, drho_darhoA, _ = computeDensity(vf, arhoA, A)
    drho_daA1 = drho_dvf * dvf_daA1

    v, dv_drho = computeSpecificVolume(rho)
    dv_daA1 = dv_drho * drho_daA1
    dv_darhoA = dv_drho * drho_darhoA

    u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)

    p = pJ
    dp_dpJ = 1

    e, de_dv, de_dp = self.eos.e(v, p)
    de_daA1 = de_dv * dv_daA1
    de_darhoA = de_dv * dv_darhoA
    de_dpJ = de_dp * dp_dpJ

    h, dh_de, dh_dp, dh_drho = computeSpecificEnthalpy(e, p, rho)
    dh_daA1 = dh_de * de_daA1 + dh_drho * drho_daA1
    dh_darhoA = dh_de * de_darhoA + dh_drho * drho_darhoA
    dh_dpJ = dh_de * de_dpJ + dh_dp * dp_dpJ

    H, dH_dh, dH_du = addKineticEnergy(h, u)
    dH_daA1 = dH_dh * dh_daA1
    dH_darhoA = dH_dh * dh_darhoA + dH_du * du_darhoA
    dH_darhouA = dH_du * du_darhouA
    dH_dpJ = dH_dh * dh_dpJ

    s, ds_dv, ds_de = self.eos.s(v, e)
    ds_daA1 = ds_dv * dv_daA1 + ds_de * de_daA1
    ds_darhoA = ds_dv * dv_darhoA + ds_de * de_darhoA
    ds_dpJ = ds_de * de_dpJ

    self.f_mass[i] = vf * rho * u * A * nx
    self.df_mass_daA1[i] = (dvf_daA1 * rho + vf * drho_daA1) * u * A * nx
    self.df_mass_darhoA[i] = vf * (drho_darhoA * u + rho * du_darhoA) * A * nx
    self.df_mass_darhouA[i] = vf * rho * du_darhouA * A * nx
    self.df_mass_dHJ[i] = 0
    self.df_mass_dsJ[i] = 0

    f_momentum = vf * (rho * u**2 + p) * A * nx
    df_momentum_daA1 = (dvf_daA1 * (rho * u**2 + p) + vf * drho_daA1 * u**2) * A * nx
    df_momentum_darhoA = vf * (drho_darhoA * u**2 + rho * 2 * u * du_darhoA) * A * nx
    df_momentum_darhouA = vf * rho * 2 * u * du_darhouA * A * nx
    df_momentum_dpJ = vf * dp_dpJ * A * nx

    self.f_energy[i] = vf * rho * H * u * A * nx
    self.df_energy_daA1[i] = (dvf_daA1 * rho * H + vf * drho_daA1 * H + vf * rho * dH_daA1) * u * A * nx
    self.df_energy_darhoA[i] = vf * (drho_darhoA * u * H + rho * du_darhoA * H + rho * u * dH_darhoA) * A * nx
    self.df_energy_darhouA[i] = vf * rho * (du_darhouA * H + u * dH_darhouA) * A * nx
    self.df_energy_dpJ[i] = vf * rho * dH_dpJ * u * A * nx
    self.df_energy_dHJ[i] = 0
    self.df_energy_dsJ[i] = 0

    self.f_entropy[i] = vf * rho * s * u * A * nx
    self.df_entropy_daA1[i] = (dvf_daA1 * rho * s + vf * drho_daA1 * s + vf * rho * ds_daA1) * u * A * nx
    self.df_entropy_darhoA[i] = vf * (drho_darhoA * u * s + rho * du_darhoA * s + rho * u * ds_darhoA) * A * nx
    self.df_entropy_darhouA[i] = vf * rho * s * du_darhouA * A * nx
    self.df_entropy_dpJ[i] = vf * rho * ds_dpJ * u * A * nx
    self.df_entropy_dHJ[i] = 0
    self.df_entropy_dsJ[i] = 0

    # add the fluxes
    i_mass = self.i_arhoA[i]
    i_momentum = self.i_arhouA[i]
    i_energy = self.i_arhoEA[i]
    i_pJ = self.i_pJ

    r[i_mass]             += self.f_mass[i]
    J[i_mass][i_mass]     += self.df_mass_darhoA[i]
    J[i_mass][i_momentum] += self.df_mass_darhouA[i]

    r[i_momentum]             += f_momentum
    J[i_momentum][i_mass]     += df_momentum_darhoA
    J[i_momentum][i_momentum] += df_momentum_darhouA
    J[i_momentum][i_pJ]       += df_momentum_dpJ

    r[i_energy]             += self.f_energy[i]
    J[i_energy][i_mass]     += self.df_energy_darhoA[i]
    J[i_energy][i_momentum] += self.df_energy_darhouA[i]
    J[i_energy][i_pJ]       += self.df_energy_dpJ[i]

  def addJunctionOutletFlux(self, i, U, r, J):
    k = self.node_indices[i]

    A = self.A[i]
    nx = self.nx[i]
    aA1 = self.dof_handler.aA1(U, k)
    arhoA = U[self.i_arhoA[i]]
    arhouA = U[self.i_arhouA[i]]
    HJ = U[self.i_HJ]
    sJ = U[self.i_sJ]

    vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

    u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)

    h, dh_dHJ, dh_du = subtractKineticEnergy(HJ, u)
    dh_darhoA  = dh_du * du_darhoA
    dh_darhouA = dh_du * du_darhouA

    p, dp_dh, dp_dsJ = self.eos.p_from_h_s(h, sJ)
    dp_darhoA  = dp_dh * dh_darhoA
    dp_darhouA = dp_dh * dh_darhouA
    dp_dHJ = dp_dh * dh_dHJ

    rho, drho_dp, drho_dsJ_partial = self.eos.rho_from_p_s(p, sJ)
    drho_darhoA  = drho_dp * dp_darhoA
    drho_darhouA = drho_dp * dp_darhouA
    drho_dHJ = drho_dp * dp_dHJ
    drho_dsJ = drho_dsJ_partial + drho_dp * dp_dsJ

    self.f_mass[i] = vf * rho * u * A * nx
    self.df_mass_daA1[i] = dvf_daA1 * rho * u * A * nx
    self.df_mass_darhoA[i]  = vf * (drho_darhoA  * u + rho * du_darhoA)  * A * nx
    self.df_mass_darhouA[i] = vf * (drho_darhouA * u + rho * du_darhouA) * A * nx
    self.df_mass_dHJ[i] = vf * drho_dHJ * u * A * nx
    self.df_mass_dsJ[i] = vf * drho_dsJ * u * A * nx

    f_momentum = vf * (rho * u**2 + p) * A * nx
    df_momentum_daA1 = dvf_daA1 * (rho * u**2 + p) * A * nx
    df_momentum_darhoA  = vf * (drho_darhoA  * u**2 + rho * 2 * u * du_darhoA  + dp_darhoA)  * A * nx
    df_momentum_darhouA = vf * (drho_darhouA * u**2 + rho * 2 * u * du_darhouA + dp_darhouA) * A * nx
    df_momentum_dHJ = vf * (drho_dHJ * u**2 + dp_dHJ) * A * nx
    df_momentum_dsJ = vf * (drho_dsJ * u**2 + dp_dsJ) * A * nx

    self.f_energy[i] = vf * rho * HJ * u * A * nx
    self.df_energy_daA1[i] = dvf_daA1 * rho * HJ * u * A * nx
    self.df_energy_darhoA[i]  = vf * (drho_darhoA  * u + rho * du_darhoA)  * HJ * A * nx
    self.df_energy_darhouA[i] = vf * (drho_darhouA * u + rho * du_darhouA) * HJ * A * nx
    self.df_energy_dpJ[i] = 0
    self.df_energy_dHJ[i] = vf * (drho_dHJ * HJ + rho) * u * A * nx
    self.df_energy_dsJ[i] = vf * drho_dsJ * HJ * u * A * nx

    self.f_entropy[i] = vf * rho * sJ * u * A * nx
    self.df_entropy_daA1[i] = dvf_daA1 * rho * sJ * u * A * nx
    self.df_entropy_darhoA[i]  = vf * (drho_darhoA  * u + rho * du_darhoA)  * sJ * A * nx
    self.df_entropy_darhouA[i] = vf * (drho_darhouA * u + rho * du_darhouA) * sJ * A * nx
    self.df_entropy_dpJ[i] = 0
    self.df_entropy_dHJ[i] = vf * drho_dHJ * sJ * u * A * nx
    self.df_entropy_dsJ[i] = vf * (drho_dsJ * sJ + rho) * u * A * nx

    # add the fluxes
    i_mass = self.i_arhoA[i]
    i_momentum = self.i_arhouA[i]
    i_energy = self.i_arhoEA[i]
    i_HJ = self.i_HJ
    i_sJ = self.i_sJ

    r[i_mass]             += self.f_mass[i]
    J[i_mass][i_mass]     += self.df_mass_darhoA[i]
    J[i_mass][i_momentum] += self.df_mass_darhouA[i]
    J[i_mass][i_HJ]       += self.df_mass_dHJ[i]
    J[i_mass][i_sJ]       += self.df_mass_dsJ[i]

    r[i_momentum]             += f_momentum
    J[i_momentum][i_mass]     += df_momentum_darhoA
    J[i_momentum][i_momentum] += df_momentum_darhouA
    J[i_momentum][i_HJ]       += df_momentum_dHJ
    J[i_momentum][i_sJ]       += df_momentum_dsJ

    r[i_energy]             += self.f_energy[i]
    J[i_energy][i_mass]     += self.df_energy_darhoA[i]
    J[i_energy][i_momentum] += self.df_energy_darhouA[i]
    J[i_energy][i_HJ]       += self.df_energy_dHJ[i]
    J[i_energy][i_sJ]       += self.df_energy_dsJ[i]

  def computeAverageQuantities(self, U):
    self.H_avg = 0
    self.s_avg = 0
    for i in xrange(self.n_meshes):
      k = self.node_indices[i]

      A = self.A[i]
      aA1 = self.dof_handler.aA1(U, k)
      arhoA = U[self.i_arhoA[i]]
      arhouA = U[self.i_arhouA[i]]
      arhoEA = U[self.i_arhoEA[i]]

      vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

      rho, drho_dvf, drho_darhoA, _ = computeDensity(vf, arhoA, A)
      drho_daA1 = drho_dvf * dvf_daA1

      u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)

      E, dE_darhoA, dE_darhoEA = computeSpecificTotalEnergy(arhoA, arhoEA)

      v, dv_drho = computeSpecificVolume(rho)
      dv_daA1 = dv_drho * drho_daA1
      dv_darhoA = dv_drho * drho_darhoA

      e, de_dE, de_du = subtractKineticEnergy(E, u)
      de_darhoA = de_dE * dE_darhoA + de_du * du_darhoA
      de_darhouA = de_du * du_darhouA
      de_darhoEA = de_dE * dE_darhoEA

      p, dp_dv, dp_de = self.eos.p(v, e)
      dp_daA1 = dp_dv * dv_daA1
      dp_darhoA = dp_dv * dv_darhoA + dp_de * de_darhoA
      dp_darhouA = dp_de * de_darhouA
      dp_darhoEA = dp_de * de_darhoEA

      s, ds_dv, ds_de = self.eos.s(v, e)
      ds_daA1 = ds_dv * dv_daA1
      ds_darhoA = ds_dv * dv_darhoA + ds_de * de_darhoA
      ds_darhouA = ds_de * de_darhouA
      ds_darhoEA = ds_de * de_darhoEA

      H, dH_dE, dH_dp, dH_drho = computeSpecificEnthalpy(E, p, rho)
      dH_daA1 = dH_dp * dp_daA1 + dH_drho * drho_daA1
      dH_darhoA = dH_dE * dE_darhoA + dH_dp * dp_darhoA + dH_drho * drho_darhoA
      dH_darhouA = dH_dp * dp_darhouA
      dH_darhoEA = dH_dE * dE_darhoEA + dH_dp * dp_darhoEA

      self.H_avg += H
      self.dH_avg_darhoA[i]  = dH_darhoA  / self.n_meshes
      self.dH_avg_darhouA[i] = dH_darhouA / self.n_meshes
      self.dH_avg_darhoEA[i] = dH_darhoEA / self.n_meshes

      self.s_avg += s
      self.ds_avg_darhoA[i]  = ds_darhoA  / self.n_meshes
      self.ds_avg_darhouA[i] = ds_darhouA / self.n_meshes
      self.ds_avg_darhoEA[i] = ds_darhoEA / self.n_meshes

    self.H_avg /= self.n_meshes
    self.s_avg /= self.n_meshes

  def applyWeaklyToNonlinearSystem(self, U_new, U_old, r, J):
    # add the boundary fluxes; characteristic theory is used to decide how
    # many pieces of information are supplied to inlets vs. outlets
    self.total_mass_flux = 0
    for i in xrange(self.n_meshes):
      A = self.A[i]
      aA1 = self.dof_handler.aA1(U_old, self.node_indices[i])
      arhoA = U_old[self.i_arhoA[i]]
      arhouA = U_old[self.i_arhouA[i]]

      u = arhouA / arhoA
      vf, _ = computeVolumeFraction(aA1, A, self.phase, self.model_type)
      rho, _, _, _ = computeDensity(vf, arhoA, A)

      self.total_mass_flux += abs(vf * rho * u * A)

      if u * self.nx[i] >= 0:
        self.addJunctionInletFlux(i, U_new, r, J)
      else:
        self.addJunctionOutletFlux(i, U_new, r, J)

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    HJ = U[self.i_HJ]
    sJ = U[self.i_sJ]

    r[self.i_constraint_mass] = sum(self.f_mass)
    J[self.i_constraint_mass,:] = 0
    for i in xrange(self.n_meshes):
      J[self.i_constraint_mass,self.i_arhoA[i]]  = self.df_mass_darhoA[i]
      J[self.i_constraint_mass,self.i_arhouA[i]] = self.df_mass_darhouA[i]
      J[self.i_constraint_mass,self.i_HJ]       += self.df_mass_dHJ[i]
      J[self.i_constraint_mass,self.i_sJ]       += self.df_mass_dsJ[i]

    if abs(self.total_mass_flux) < 1e-12:
      # replace energy and entropy balance constraints with averages
      self.computeAverageQuantities(U)
      r[self.i_constraint_energy]  = HJ - self.H_avg
      r[self.i_constraint_entropy] = sJ - self.s_avg
      J[self.i_constraint_energy,:]  = 0
      J[self.i_constraint_entropy,:] = 0
      for i in xrange(self.n_meshes):
        J[self.i_constraint_energy,self.i_arhoA[i]]  = -self.dH_avg_darhoA[i]
        J[self.i_constraint_energy,self.i_arhouA[i]] = -self.dH_avg_darhouA[i]
        J[self.i_constraint_energy,self.i_arhoEA[i]] = -self.dH_avg_darhoEA[i]
        J[self.i_constraint_energy,self.i_HJ] = 1

        J[self.i_constraint_entropy,self.i_arhoA[i]]  = -self.ds_avg_darhoA[i]
        J[self.i_constraint_entropy,self.i_arhouA[i]] = -self.ds_avg_darhouA[i]
        J[self.i_constraint_entropy,self.i_arhoEA[i]] = -self.ds_avg_darhoEA[i]
        J[self.i_constraint_entropy,self.i_sJ] = 1
    else:
      r[self.i_constraint_energy]  = sum(self.f_energy)
      r[self.i_constraint_entropy] = sum(self.f_entropy)
      J[self.i_constraint_energy,:]  = 0
      J[self.i_constraint_entropy,:] = 0
      for i in xrange(self.n_meshes):
        J[self.i_constraint_energy,self.i_arhoA[i]]  = self.df_energy_darhoA[i]
        J[self.i_constraint_energy,self.i_arhouA[i]] = self.df_energy_darhouA[i]
        J[self.i_constraint_energy,self.i_pJ]       += self.df_energy_dpJ[i]
        J[self.i_constraint_energy,self.i_HJ]       += self.df_energy_dHJ[i]
        J[self.i_constraint_energy,self.i_sJ]       += self.df_energy_dsJ[i]

        J[self.i_constraint_entropy,self.i_arhoA[i]]  = self.df_entropy_darhoA[i]
        J[self.i_constraint_entropy,self.i_arhouA[i]] = self.df_entropy_darhouA[i]
        J[self.i_constraint_entropy,self.i_pJ]       += self.df_entropy_dpJ[i]
        J[self.i_constraint_entropy,self.i_HJ]       += self.df_entropy_dHJ[i]
        J[self.i_constraint_entropy,self.i_sJ]       += self.df_entropy_dsJ[i]

  def applyStronglyToLinearSystemMatrix(self, A):
    pass

  def applyStronglyToLinearSystemRHSVector(self, U_old, b):
    pass
