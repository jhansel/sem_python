from Junction import Junction, JunctionParameters
from enums import ModelType, VariableName
from thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy, computeSpecificEnthalpy

class Junction1PhaseParameters(JunctionParameters):
  def __init__(self):
    JunctionParameters.__init__(self)
    self.registerIntParameter("phase", "Index of phase to which BC is applied")

## Base class for 1-phase junctions
class Junction1Phase(Junction):
  def __init__(self, params):
    Junction.__init__(self, params)
    self.phase = params.get("phase")
    self.eos = self.eos_list[self.phase]

    # get variable indices
    if (self.model_type == ModelType.TwoPhase):
      self.aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]
    self.arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][self.phase]
    self.arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][self.phase]
    self.arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][self.phase]

    # initialize lists
    self.aA1 = [0] * self.n_meshes
    self.arhoA = [0] * self.n_meshes
    self.arhouA = [0] * self.n_meshes
    self.arhoEA = [0] * self.n_meshes

    self.vf = [0] * self.n_meshes
    self.dvf_daA1 = [0] * self.n_meshes

    self.rho = [0] * self.n_meshes
    self.drho_daA1 = [0] * self.n_meshes
    self.drho_darhoA = [0] * self.n_meshes

    self.v = [0] * self.n_meshes
    self.dv_daA1 = [0] * self.n_meshes
    self.dv_darhoA = [0] * self.n_meshes

    self.u = [0] * self.n_meshes
    self.du_darhoA = [0] * self.n_meshes
    self.du_darhouA = [0] * self.n_meshes

    self.E = [0] * self.n_meshes
    self.dE_darhoA = [0] * self.n_meshes
    self.dE_darhoEA = [0] * self.n_meshes

    self.e = [0] * self.n_meshes
    self.de_darhoA = [0] * self.n_meshes
    self.de_darhouA = [0] * self.n_meshes
    self.de_darhoEA = [0] * self.n_meshes

    self.p = [0] * self.n_meshes
    self.dp_daA1 = [0] * self.n_meshes
    self.dp_darhoA = [0] * self.n_meshes
    self.dp_darhouA = [0] * self.n_meshes
    self.dp_darhoEA = [0] * self.n_meshes

    self.c = [0] * self.n_meshes
    self.dc_daA1 = [0] * self.n_meshes
    self.dc_darhoA = [0] * self.n_meshes
    self.dc_darhouA = [0] * self.n_meshes
    self.dc_darhoEA = [0] * self.n_meshes

    self.s = [0] * self.n_meshes
    self.ds_daA1 = [0] * self.n_meshes
    self.ds_darhoA = [0] * self.n_meshes
    self.ds_darhouA = [0] * self.n_meshes
    self.ds_darhoEA = [0] * self.n_meshes

    self.h = [0] * self.n_meshes
    self.dh_daA1 = [0] * self.n_meshes
    self.dh_darhoA = [0] * self.n_meshes
    self.dh_darhouA = [0] * self.n_meshes
    self.dh_darhoEA = [0] * self.n_meshes

    self.h0 = [0] * self.n_meshes
    self.dh0_daA1 = [0] * self.n_meshes
    self.dh0_darhoA = [0] * self.n_meshes
    self.dh0_darhouA = [0] * self.n_meshes
    self.dh0_darhoEA = [0] * self.n_meshes

    self.p0 = [0] * self.n_meshes
    self.dp0_daA1 = [0] * self.n_meshes
    self.dp0_darhoA = [0] * self.n_meshes
    self.dp0_darhouA = [0] * self.n_meshes
    self.dp0_darhoEA = [0] * self.n_meshes

    self.f_mass = [0] * self.n_meshes
    self.df_mass_darhouA = [0] * self.n_meshes

    self.f_momentum = [0] * self.n_meshes
    self.df_momentum_daA1 = [0] * self.n_meshes
    self.df_momentum_darhoA = [0] * self.n_meshes
    self.df_momentum_darhouA = [0] * self.n_meshes
    self.df_momentum_darhoEA = [0] * self.n_meshes

    self.f_energy = [0] * self.n_meshes
    self.df_energy_daA1 = [0] * self.n_meshes
    self.df_energy_darhoA = [0] * self.n_meshes
    self.df_energy_darhouA = [0] * self.n_meshes
    self.df_energy_darhoEA = [0] * self.n_meshes

  def setDoFIndices(self):
    # indices for all meshes
    self.i_arhoA = [self.dof_handler.i(k, self.arhoA_index) for k in self.node_indices]
    self.i_arhouA = [self.dof_handler.i(k, self.arhouA_index) for k in self.node_indices]
    self.i_arhoEA = [self.dof_handler.i(k, self.arhoEA_index) for k in self.node_indices]

  def computeFlowQuantities(self, U):
    for i in xrange(self.n_meshes):
      k = self.node_indices[i]

      aA1 = self.dof_handler.aA1(U, k)
      arhoA = U[self.i_arhoA[i]]
      arhouA = U[self.i_arhouA[i]]
      arhoEA = U[self.i_arhoEA[i]]

      self.vf[i], self.dvf_daA1[i] = computeVolumeFraction(aA1, self.A[i], self.phase, self.model_type)

      self.rho[i], drho_dvf, self.drho_darhoA[i], _ = computeDensity(self.vf[i], arhoA, self.A[i])
      self.drho_daA1[i] = drho_dvf * self.dvf_daA1[i]

      self.v[i], dv_drho = computeSpecificVolume(self.rho[i])
      self.dv_daA1[i] = dv_drho * self.drho_daA1[i]
      self.dv_darhoA[i] = dv_drho * self.drho_darhoA[i]

      self.u[i], self.du_darhoA[i], self.du_darhouA[i] = computeVelocity(arhoA, arhouA)

      self.E[i], self.dE_darhoA[i], self.dE_darhoEA[i] = computeSpecificTotalEnergy(arhoA, arhoEA)

      self.e[i], de_du, de_dE = computeSpecificInternalEnergy(self.u[i], self.E[i])
      self.de_darhoA[i] = de_du * self.du_darhoA[i] + de_dE * self.dE_darhoA[i]
      self.de_darhouA[i] = de_du * self.du_darhouA[i]
      self.de_darhoEA[i] = de_dE * self.dE_darhoEA[i]

      self.p[i], dp_dv, dp_de = self.eos.p(self.v[i], self.e[i])
      self.dp_daA1[i] = dp_dv * self.dv_daA1[i]
      self.dp_darhoA[i] = dp_dv * self.dv_darhoA[i] + dp_de * self.de_darhoA[i]
      self.dp_darhouA[i] = dp_de * self.de_darhouA[i]
      self.dp_darhoEA[i] = dp_de * self.de_darhoEA[i]

      self.c[i], dc_dv, dc_de = self.eos.c(self.v[i], self.e[i])
      self.dc_daA1[i] = dc_dv * self.dv_daA1[i]
      self.dc_darhoA[i] = dc_dv * self.dv_darhoA[i] + dc_de * self.de_darhoA[i]
      self.dc_darhouA[i] = dc_de * self.de_darhouA[i]
      self.dc_darhoEA[i] = dc_de * self.de_darhoEA[i]

      self.s[i], ds_dv, ds_de = self.eos.s(self.v[i], self.e[i])
      self.ds_daA1[i] = ds_dv * self.dv_daA1[i]
      self.ds_darhoA[i] = ds_dv * self.dv_darhoA[i] + ds_de * self.de_darhoA[i]
      self.ds_darhouA[i] = ds_de * self.de_darhouA[i]
      self.ds_darhoEA[i] = ds_de * self.de_darhoEA[i]

      self.h[i], dh_de, dh_dp, dh_drho = computeSpecificEnthalpy(self.e[i], self.p[i], self.rho[i])
      self.dh_daA1[i] = dh_dp * self.dp_daA1[i]
      self.dh_darhoA[i] = dh_de * self.de_darhoA[i] + dh_dp * self.dp_darhoA[i] + dh_drho * self.drho_darhoA[i]
      self.dh_darhouA[i] = dh_de * self.de_darhouA[i] + dh_dp * self.dp_darhouA[i]
      self.dh_darhoEA[i] = dh_de * self.de_darhoEA[i] + dh_dp * self.dp_darhoEA[i]

      self.h0[i] = self.h[i] + 0.5 * self.u[i]**2
      self.dh0_daA1[i] = self.dh_daA1[i]
      self.dh0_darhoA[i] = self.dh_darhoA[i] + self.u[i] * self.du_darhoA[i]
      self.dh0_darhouA[i] = self.dh_darhouA[i] + self.u[i] * self.du_darhouA[i]
      self.dh0_darhoEA[i] = self.dh_darhoEA[i]

      self.p0[i], dp0_dh0, dp0_ds = self.eos.p_from_h_s(self.h0[i], self.s[i])
      self.dp0_daA1[i] = dp0_dh0 * self.dh0_daA1[i] + dp0_ds * self.ds_daA1[i]
      self.dp0_darhoA[i] = dp0_dh0 * self.dh0_darhoA[i] + dp0_ds * self.ds_darhoA[i]
      self.dp0_darhouA[i] = dp0_dh0 * self.dh0_darhouA[i] + dp0_ds * self.ds_darhouA[i]
      self.dp0_darhoEA[i] = dp0_dh0 * self.dh0_darhoEA[i] + dp0_ds * self.ds_darhoEA[i]

  def computeFluxes(self, U):
    for i in xrange(self.n_meshes):
      arhouA = U[self.i_arhouA[i]]
      arhoEA = U[self.i_arhoEA[i]]

      nx = self.nx[i]

      vf = self.vf[i]
      u = self.u[i]
      p = self.p[i]

      self.f_mass[i] = arhouA * nx
      self.df_mass_darhouA[i] = nx

      self.f_momentum[i] = (arhouA * u + vf * p) * nx
      self.df_momentum_daA1[i] = (self.dvf_daA1[i] * p + vf * self.dp_daA1[i]) * nx
      self.df_momentum_darhoA[i] = (arhouA * self.du_darhoA[i] + vf * self.dp_darhoA[i]) * nx
      self.df_momentum_darhouA[i] = (arhouA * self.du_darhouA[i] + u + vf * self.dp_darhouA[i]) * nx
      self.df_momentum_darhoEA[i] = vf * self.dp_darhoEA[i] * nx

      self.f_energy[i] = (arhoEA + vf * p) * u * nx
      self.df_energy_daA1[i] = (self.dvf_daA1[i] * p + vf * self.dp_daA1[i]) * u * nx
      self.df_energy_darhoA[i] = ((arhoEA + vf * p) * self.du_darhoA[i] + vf * self.dp_darhoA[i] * u) * nx
      self.df_energy_darhouA[i] = ((arhoEA + vf * p) * self.du_darhouA[i] + vf * self.dp_darhouA[i] * u) * nx
      self.df_energy_darhoEA[i] = (1 + vf * self.dp_darhoEA[i]) * u * nx

  def applyWeaklyToNonlinearSystem(self, U, U_old, r, J):
    self.computeFlowQuantities(U)
    self.computeFluxes(U)
    for m in xrange(self.n_meshes):
      r[self.i_arhoA[m]] += self.f_mass[m]
      r[self.i_arhouA[m]] += self.f_momentum[m]
      r[self.i_arhoEA[m]] += self.f_energy[m]

      J[self.i_arhoA[m],self.i_arhouA[m]] += self.df_mass_darhouA[m]
      J[self.i_arhouA[m],self.i_arhoA[m]] += self.df_momentum_darhoA[m]
      J[self.i_arhouA[m],self.i_arhouA[m]] += self.df_momentum_darhouA[m]
      J[self.i_arhouA[m],self.i_arhoEA[m]] += self.df_momentum_darhoEA[m]
      J[self.i_arhoEA[m],self.i_arhoA[m]] += self.df_energy_darhoA[m]
      J[self.i_arhoEA[m],self.i_arhouA[m]] += self.df_energy_darhouA[m]
      J[self.i_arhoEA[m],self.i_arhoEA[m]] += self.df_energy_darhoEA[m]

      if self.model_type == ModelType.TwoPhase:
        J[self.i_arhouA[m],self.i_aA1[m]] += self.df_momentum_daA1[m]
        J[self.i_arhoEA[m],self.i_aA1[m]] += self.df_energy_daA1[m]
