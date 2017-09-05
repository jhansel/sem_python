from Junction1Phase import Junction1Phase, Junction1PhaseParameters
from thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy
from enums import ModelType
from error_utilities import error

class FluxJunctionParameters(Junction1PhaseParameters):
  def __init__(self):
    Junction1PhaseParameters.__init__(self)

## Base class for junctions that compute boundary fluxes
class FluxJunction(Junction1Phase):
  def __init__(self, params):
    Junction1Phase.__init__(self, params)
    if self.n_meshes != 2:
      error("FluxJunction is only implemented for connecting 2 meshes.")

  def setDoFIndices(self):
    Junction1Phase.setDoFIndices(self)

  def computeFluxes(self, U):
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

    for i in xrange(self.n_meshes):
      nx = self.nx[i]

      A = self.dof_handler.A[self.node_indices[i]]
      aA1 = self.dof_handler.aA1(U, self.node_indices[i])
      vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)
      arhoA = U[self.i_arhoA[i]]
      arhouA = U[self.i_arhouA[i]]
      arhoEA = U[self.i_arhoEA[i]]

      u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)

      rho, drho_dvf, drho_darhoA, _ = computeDensity(vf, arhoA, A)
      drho_daA1 = drho_dvf * dvf_daA1

      v, dv_drho = computeSpecificVolume(rho)
      dv_daA1 = dv_drho * drho_daA1
      dv_darhoA = dv_drho * drho_darhoA

      E, dE_darhoA, dE_darhoEA = computeSpecificTotalEnergy(arhoA, arhoEA)

      e, de_du, de_dE = computeSpecificInternalEnergy(u, E)
      de_darhoA = de_du * du_darhoA + de_dE * dE_darhoA
      de_darhouA = de_du * du_darhouA
      de_darhoEA = de_dE * dE_darhoEA

      p, dp_dv, dp_de = self.eos.p(v, e)
      dp_daA1 = dp_dv * dv_daA1
      dp_darhoA = dp_dv * dv_darhoA + dp_de * de_darhoA
      dp_darhouA = dp_de * de_darhouA
      dp_darhoEA = dp_de * de_darhoEA

      self.f_mass[i] = arhouA * nx
      self.df_mass_darhouA[i] = nx

      self.f_momentum[i] = (arhouA * u + vf * p) * nx
      self.df_momentum_daA1[i] = (dvf_daA1 * p + vf * dp_daA1) * nx
      self.df_momentum_darhoA[i] = (arhouA * du_darhoA + vf * dp_darhoA) * nx
      self.df_momentum_darhouA[i] = (arhouA * du_darhouA + u + vf * dp_darhouA) * nx
      self.df_momentum_darhoEA[i] = vf * dp_darhoEA * nx

      self.f_energy[i] = (arhoEA + vf * p) * u * nx
      self.df_energy_daA1[i] = (dvf_daA1 * p + vf * dp_daA1) * u * nx
      self.df_energy_darhoA[i] = ((arhoEA + vf * p) * du_darhoA + vf * dp_darhoA * u) * nx
      self.df_energy_darhouA[i] = ((arhoEA + vf * p) * du_darhouA + vf * dp_darhouA * u) * nx
      self.df_energy_darhoEA[i] = (1 + vf * dp_darhoEA) * u * nx
