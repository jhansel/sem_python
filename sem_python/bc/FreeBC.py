from ..base.enums import ModelType
from .OnePhaseBC import OnePhaseBC, OnePhaseBCParameters
from ..closures.thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeSpecificVolume, computeVelocity, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy

class FreeBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)

class FreeBC(OnePhaseBC):
  def __init__(self, params):
    OnePhaseBC.__init__(self, params)

  def applyWeakBC(self, U, r, J):
    A = self.dof_handler.A[self.k]
    aA1 = self.dof_handler.aA1(U, self.k)
    arhoA = U[self.i_arhoA]
    arhouA = U[self.i_arhouA]
    arhoEA = U[self.i_arhoEA]

    vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

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

    # mass
    r[self.i_arhoA] += arhouA * self.nx
    J[self.i_arhoA,self.i_arhouA] += self.nx

    # momentum
    r[self.i_arhouA] += (arhouA * u + vf * p * A) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhouA,self.i_aA1] += (dvf_daA1 * p + vf * dp_daA1) * A * self.nx
    J[self.i_arhouA,self.i_arhoA] += (arhouA * du_darhoA + vf * dp_darhoA * A) * self.nx
    J[self.i_arhouA,self.i_arhouA] += (arhouA * du_darhouA + u + vf * dp_darhouA * A) * self.nx
    J[self.i_arhouA,self.i_arhoEA] += vf * dp_darhoEA * A * self.nx

    # energy
    r[self.i_arhoEA] += (arhoEA + vf * p * A) * u * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoEA,self.i_aA1] += (dvf_daA1 * p + vf * dp_daA1) * A * u * self.nx
    J[self.i_arhoEA,self.i_arhoA] += ((arhoEA + vf * p * A) * du_darhoA + vf * dp_darhoA * A * u) * self.nx
    J[self.i_arhoEA,self.i_arhouA] += ((arhoEA + vf * p * A) * du_darhouA + vf * dp_darhouA * A * u) * self.nx
    J[self.i_arhoEA,self.i_arhoEA] += (1 + vf * dp_darhoEA * A) * u * self.nx
