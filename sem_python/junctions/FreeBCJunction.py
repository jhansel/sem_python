from .FluxJunction import FluxJunction, FluxJunctionParameters
from ..base.enums import ModelType

class FreeBCJunctionParameters(FluxJunctionParameters):
  def __init__(self):
    FluxJunctionParameters.__init__(self)

## Junction that adds boundary fluxes to each mesh, computed from their own solutions
class FreeBCJunction(FluxJunction):
  def __init__(self, params):
    FluxJunction.__init__(self, params)

  def setDoFIndices(self):
    FluxJunction.setDoFIndices(self)

  def applyWeaklyToNonlinearSystem(self, U, U_old, r, J):
    self.computeFluxes(U)
    for m in xrange(self.n_meshes):
      r[self.i_arho[m]] += self.f_mass[m]
      r[self.i_arhou[m]] += self.f_momentum[m]
      r[self.i_arhoE[m]] += self.f_energy[m]

      J[self.i_arho[m],self.i_arhou[m]] += self.df_mass_darhou[m]
      J[self.i_arhou[m],self.i_arho[m]] += self.df_momentum_darho[m]
      J[self.i_arhou[m],self.i_arhou[m]] += self.df_momentum_darhou[m]
      J[self.i_arhou[m],self.i_arhoE[m]] += self.df_momentum_darhoE[m]
      J[self.i_arhoE[m],self.i_arho[m]] += self.df_energy_darho[m]
      J[self.i_arhoE[m],self.i_arhou[m]] += self.df_energy_darhou[m]
      J[self.i_arhoE[m],self.i_arhoE[m]] += self.df_energy_darhoE[m]

      if self.model_type == ModelType.TwoPhase:
        J[self.i_arhou[m],self.i_vf1[m]] += self.df_momentum_dvf1[m]
        J[self.i_arhoE[m],self.i_vf1[m]] += self.df_energy_dvf1[m]

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    pass

  def applyStronglyToLinearSystemMatrix(self, A):
    error("Not implemented")

  def applyStronglyToLinearSystemRHSVector(self, U_old, b):
    error("Not implemented")
