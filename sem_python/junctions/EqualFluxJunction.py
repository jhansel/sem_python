from FreeBCJunction import FreeBCJunction, FreeBCJunctionParameters
from enums import ModelType

class EqualFluxJunctionParameters(FreeBCJunctionParameters):
  def __init__(self):
    FreeBCJunctionParameters.__init__(self)

## Junction that enforces that the fluxes are equal at the junction
#
#  The first mesh will have its residuals replaced strongly, and the second
#  mesh will add a boundary flux computed with its solution.
class EqualFluxJunction(FreeBCJunction):
  def __init__(self, params):
    FreeBCJunction.__init__(self, params)

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    m = 0
    r[self.i_arho[m]] = sum(self.f_mass)
    r[self.i_arhou[m]] = sum(self.f_momentum)
    r[self.i_arhoE[m]] = sum(self.f_energy)
    J[self.i_arho[m],:] = 0
    J[self.i_arhou[m],:] = 0
    J[self.i_arhoE[m],:] = 0

    for n in xrange(self.n_meshes):
      J[self.i_arho[m],self.i_arhou[n]] = self.df_mass_darhou[n]

      J[self.i_arhou[m],self.i_arho[n]] = self.df_momentum_darho[n]
      J[self.i_arhou[m],self.i_arhou[n]] = self.df_momentum_darhou[n]
      J[self.i_arhou[m],self.i_arhoE[n]] = self.df_momentum_darhoE[n]

      J[self.i_arhoE[m],self.i_arho[n]] = self.df_energy_darho[n]
      J[self.i_arhoE[m],self.i_arhou[n]] = self.df_energy_darhou[n]
      J[self.i_arhoE[m],self.i_arhoE[n]] = self.df_energy_darhoE[n]

      if self.model_type == ModelType.TwoPhase:
        J[self.i_arhou[m],self.i_vf1[n]] = self.df_momentum_dvf1[n]
        J[self.i_arhoE[m],self.i_vf1[n]] = self.df_energy_dvf1[n]
