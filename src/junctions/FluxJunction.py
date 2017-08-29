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

  def computeFluxes(self, U):
    self.f_mass = [0] * self.n_meshes
    self.df_mass_darhou = [0] * self.n_meshes

    self.f_momentum = [0] * self.n_meshes
    self.df_momentum_dvf1 = [0] * self.n_meshes
    self.df_momentum_darho = [0] * self.n_meshes
    self.df_momentum_darhou = [0] * self.n_meshes
    self.df_momentum_darhoE = [0] * self.n_meshes

    self.f_energy = [0] * self.n_meshes
    self.df_energy_dvf1 = [0] * self.n_meshes
    self.df_energy_darho = [0] * self.n_meshes
    self.df_energy_darhou = [0] * self.n_meshes
    self.df_energy_darhoE = [0] * self.n_meshes

    for i in xrange(self.n_meshes):
      nx = self.nx[i]

      vf1 = self.dof_handler.getVolumeFraction(U, self.node_indices[i])
      vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
      arho = U[self.i_arho[i]]
      arhou = U[self.i_arhou[i]]
      arhoE = U[self.i_arhoE[i]]

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

      p, dp_dv, dp_de = self.eos.p(v, e)
      dp_dvf1 = dp_dv * dv_dvf1
      dp_darho = dp_dv * dv_darho + dp_de * de_darho
      dp_darhou = dp_de * de_darhou
      dp_darhoE = dp_de * de_darhoE

      self.f_mass[i] = arhou * nx
      self.df_mass_darhou[i] = nx

      self.f_momentum[i] = (arhou * u + vf * p) * nx
      self.df_momentum_dvf1[i] = (dvf_dvf1 * p + vf * dp_dvf1) * nx
      self.df_momentum_darho[i] = (arhou * du_darho + vf * dp_darho) * nx
      self.df_momentum_darhou[i] = (arhou * du_darhou + u + vf * dp_darhou) * nx
      self.df_momentum_darhoE[i] = vf * dp_darhoE * nx

      self.f_energy[i] = (arhoE + vf * p) * u * nx
      self.df_energy_dvf1[i] = (dvf_dvf1 * p + vf * dp_dvf1) * u * nx
      self.df_energy_darho[i] = ((arhoE + vf * p) * du_darho + vf * dp_darho * u) * nx
      self.df_energy_darhou[i] = ((arhoE + vf * p) * du_darhou + vf * dp_darhou * u) * nx
      self.df_energy_darhoE[i] = (1 + vf * dp_darhoE) * u * nx
