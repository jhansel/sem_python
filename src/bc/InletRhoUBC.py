from enums import ModelType
from OnePhaseBC import OnePhaseBC, OnePhaseBCParameters
from thermodynamic_functions import computeVolumeFraction, computeSpecificVolume, \
  computeSpecificTotalEnergy, computeSpecificInternalEnergy

class InletRhoUBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)
    self.registerFloatParameter("rho", "Specified density")
    self.registerFloatParameter("u", "Specified velocity")

class InletRhoUBC(OnePhaseBC):
  def __init__(self, params):
    OnePhaseBC.__init__(self, params)
    self.rho = params.get("rho")
    self.u = params.get("u")

  def applyWeakBC(self, U, r, J):
    vf1 = self.dof_handler.getVolumeFraction(U, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arhoE = U[self.i_arhoE]

    arhoBC = vf * self.rho
    darhoBC_dvf1 = self.rho * dvf_dvf1

    arhouBC = vf * self.rho * self.u
    darhouBC_dvf1 = self.rho * self.u * dvf_dvf1

    v, _ = computeSpecificVolume(self.rho)

    E, dE_darhoBC, dE_darhoE = computeSpecificTotalEnergy(arhoBC, arhoE)
    dE_dvf1 = dE_darhoBC * darhoBC_dvf1

    e, _, de_dE = computeSpecificInternalEnergy(self.u, E)
    de_dvf1 = de_dE * dE_dvf1
    de_darhoE = de_dE * dE_darhoE

    p, _, dp_de = self.eos.p(v, e)
    dp_dvf1 = dp_de * de_dvf1
    dp_darhoE = dp_de * de_darhoE

    # momentum
    r[self.i_arhou] += (arhouBC * self.u + vf * p) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhou,self.i_vf1] += (self.u * darhouBC_dvf1 + dvf_dvf1 * p + vf * dp_dvf1) * self.nx
    J[self.i_arhou,self.i_arhoE] += vf * dp_darhoE * self.nx

    # energy
    r[self.i_arhoE] += self.u * (arhoE + vf * p) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoE,self.i_vf1] += self.u * (dvf_dvf1 * p + vf * dp_dvf1) * self.nx
    J[self.i_arhoE,self.i_arhoE] += self.u * (1 + vf * dp_darhoE) * self.nx

  def applyStrongBCNonlinearSystem(self, U, r, J):
    vf1 = self.dof_handler.getVolumeFraction(U, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arho = U[self.i_arho]

    arhoBC = vf * self.rho
    darhoBC_dvf1 = self.rho * dvf_dvf1

    # mass
    r[self.i_arho] = arho - arhoBC
    J[self.i_arho,:] = 0
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arho,self.i_vf1] = - darhoBC_dvf1
    J[self.i_arho,self.i_arho] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_arho,:] = 0
    A[self.i_arho,self.i_arho] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    vf1 = self.dof_handler.getVolumeFraction(U_old, self.k)
    vf, dvf_dvf1 = computeVolumeFraction(vf1, self.phase, self.model_type)
    arho = vf * self.rho
    b[self.i_arho] = arho
