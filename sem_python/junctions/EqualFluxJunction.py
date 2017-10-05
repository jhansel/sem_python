from .Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..base.enums import ModelType

class EqualFluxJunctionParameters(Junction1PhaseParameters):
  def __init__(self):
    Junction1PhaseParameters.__init__(self)

## Junction that enforces that the fluxes are equal at the junction
#
#  The first mesh will have its residuals replaced strongly, and the second
#  mesh will add a boundary flux computed with its solution.
class EqualFluxJunction(Junction1Phase):
  def __init__(self, params):
    Junction1Phase.__init__(self, params)

    self.initializeVariableVectors()

  def applyStronglyToNonlinearSystem(self, U, U_old, r, J):
    m = 0
    r[self.i_arhoA[m]] = sum(self.f_mass)
    r[self.i_arhouA[m]] = sum(self.f_momentum)
    r[self.i_arhoEA[m]] = sum(self.f_energy)
    J[self.i_arhoA[m],:] = 0
    J[self.i_arhouA[m],:] = 0
    J[self.i_arhoEA[m],:] = 0

    for n in xrange(self.n_meshes):
      J[self.i_arhoA[m],self.i_arhouA[n]] = self.df_mass_darhouA[n]

      J[self.i_arhouA[m],self.i_arhoA[n]] = self.df_momentum_darhoA[n]
      J[self.i_arhouA[m],self.i_arhouA[n]] = self.df_momentum_darhouA[n]
      J[self.i_arhouA[m],self.i_arhoEA[n]] = self.df_momentum_darhoEA[n]

      J[self.i_arhoEA[m],self.i_arhoA[n]] = self.df_energy_darhoA[n]
      J[self.i_arhoEA[m],self.i_arhouA[n]] = self.df_energy_darhouA[n]
      J[self.i_arhoEA[m],self.i_arhoEA[n]] = self.df_energy_darhoEA[n]

      if self.model_type == ModelType.TwoPhase:
        J[self.i_arhouA[m],self.i_aA1[n]] = self.df_momentum_daA1[n]
        J[self.i_arhoEA[m],self.i_aA1[n]] = self.df_energy_daA1[n]
