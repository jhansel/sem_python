from ..base.enums import ModelType
from .OnePhaseBC import OnePhaseBC, OnePhaseBCParameters
from ..closures.thermodynamic_functions import computeVolumeFraction, \
  computeDensity, computeVelocity, computeSpecificVolume

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
    arho = U[self.i_arho]
    arhou = U[self.i_arhou]

    u, du_darho, du_darhou = computeVelocity(arho, arhou)

    # mass
    r[self.i_arho] += arhou * self.nx
    J[self.i_arho,self.i_arhou] += self.nx

    # momentum
    r[self.i_arhou] += (arhou * u + vf * self.p) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhou,self.i_vf1] += dvf_dvf1 * self.p * self.nx
    J[self.i_arhou,self.i_arho] += arhou * du_darho * self.nx
    J[self.i_arhou,self.i_arhou] += (arhou * du_darhou + u) * self.nx

  def applyStrongBCNonlinearSystem(self, U, r, J):
    vf1 = self.dof_handler.getVolumeFraction(U, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arho = U[self.i_arho]
    arhou = U[self.i_arhou]
    arhoE_solution = U[self.i_arhoE]

    u, du_darho, du_darhou = computeVelocity(arho, arhou)

    rho, drho_dvf, drho_darho = computeDensity(vf, arho)
    drho_dvf1 = drho_dvf * dvf_dvf1

    v, dv_drho = computeSpecificVolume(rho)
    dv_dvf1 = dv_drho * drho_dvf1
    dv_darho = dv_drho * drho_darho

    e, de_dv, de_dp = self.eos.e(v, self.p)
    de_dvf1 = de_dv * dv_dvf1
    de_darho = de_dv * dv_darho

    arhoE = arho * (e + 0.5 * u * u)
    darhoE_dvf1 = arho * de_dv * dv_dvf1
    darhoE_darho = (e + 0.5 * u * u) + arho * (de_dv * dv_darho + u * du_darho)
    darhoE_darhou = arho * u * du_darhou

    # energy
    r[self.i_arhoE] = arhoE_solution - arhoE
    J[self.i_arhoE,:] = 0
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoE,self.i_vf1] = - darhoE_dvf1
    J[self.i_arhoE,self.i_arho] = - darhoE_darho
    J[self.i_arhoE,self.i_arhou] = - darhoE_darhou
    J[self.i_arhoE,self.i_arhoE] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_arhoE,:] = 0
    A[self.i_arhoE,self.i_arhoE] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    vf1 = self.dof_handler.getVolumeFraction(U_old, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arho = U_old[self.i_arho]
    arhou = U_old[self.i_arhou]

    u, du_darho, du_darhou = computeVelocity(arho, arhou)
    rho, drho_dvf, drho_darho = computeDensity(vf, arho)
    v, dv_drho = computeSpecificVolume(rho)
    e, de_dv, de_dp = self.eos.e(v, self.p)
    arhoE = arho * (e + 0.5 * u * u)

    b[self.i_arhoE] = arhoE
