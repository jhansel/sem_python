from FluxJunction import FluxJunction, FluxJunctionParameters
from enums import ModelType

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

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    pass

  def applyStronglyToLinearSystemMatrix(self, A):
    error("Not implemented")

  def applyStronglyToLinearSystemRHSVector(self, U_old, b):
    error("Not implemented")
