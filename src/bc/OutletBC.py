from enums import ModelType
from OnePhaseBC import OnePhaseBC, OnePhaseBCParameters
from thermodynamic_functions import computeVolumeFraction, computeDensity, \
  computeVelocity, computeSpecificVolume

## Parameters class for OutletBC
class OutletBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)
    self.registerFloatParameter("p", "specified outlet pressure")

class OutletBC(OnePhaseBC):
  def __init__(self, params):
    OnePhaseBC.__init__(self, params)
    self.p = params.get("p")

  def applyWeakBC(self, U, r, J):
    vf1 = self.dof_handler.getVolumeFraction(U, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arhoA = U[self.i_arhoA]
    arhouA = U[self.i_arhouA]

    u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)

    # mass
    r[self.i_arhoA] += arhouA * self.nx
    J[self.i_arhoA,self.i_arhouA] += self.nx

    # momentum
    r[self.i_arhouA] += (arhouA * u + vf * self.p) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhouA,self.i_vf1] += dvf_dvf1 * self.p * self.nx
    J[self.i_arhouA,self.i_arhoA] += arhouA * du_darhoA * self.nx
    J[self.i_arhouA,self.i_arhouA] += (arhouA * du_darhouA + u) * self.nx

  def applyStrongBCNonlinearSystem(self, U, r, J):
    vf1 = self.dof_handler.getVolumeFraction(U, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arhoA = U[self.i_arhoA]
    arhouA = U[self.i_arhouA]
    arhoEA_solution = U[self.i_arhoEA]

    u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)

    rho, drho_dvf, drho_darhoA = computeDensity(vf, arhoA)
    drho_dvf1 = drho_dvf * dvf_dvf1

    v, dv_drho = computeSpecificVolume(rho)
    dv_dvf1 = dv_drho * drho_dvf1
    dv_darhoA = dv_drho * drho_darhoA

    e, de_dv, de_dp = self.eos.e(v, self.p)
    de_dvf1 = de_dv * dv_dvf1
    de_darhoA = de_dv * dv_darhoA

    arhoEA = arhoA * (e + 0.5 * u * u)
    darhoEA_dvf1 = arhoA * de_dv * dv_dvf1
    darhoEA_darhoA = (e + 0.5 * u * u) + arhoA * (de_dv * dv_darhoA + u * du_darhoA)
    darhoEA_darhouA = arhoA * u * du_darhouA

    # energy
    r[self.i_arhoEA] = arhoEA_solution - arhoEA
    J[self.i_arhoEA,:] = 0
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoEA,self.i_vf1] = - darhoEA_dvf1
    J[self.i_arhoEA,self.i_arhoA] = - darhoEA_darhoA
    J[self.i_arhoEA,self.i_arhouA] = - darhoEA_darhouA
    J[self.i_arhoEA,self.i_arhoEA] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_arhoEA,:] = 0
    A[self.i_arhoEA,self.i_arhoEA] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    vf1 = self.dof_handler.getVolumeFraction(U_old, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arhoA = U_old[self.i_arhoA]
    arhouA = U_old[self.i_arhouA]

    u, du_darhoA, du_darhouA = computeVelocity(arhoA, arhouA)
    rho, drho_dvf, drho_darhoA = computeDensity(vf, arhoA)
    v, dv_drho = computeSpecificVolume(rho)
    e, de_dv, de_dp = self.eos.e(v, self.p)
    arhoEA = arhoA * (e + 0.5 * u * u)

    b[self.i_arhoEA] = arhoEA
