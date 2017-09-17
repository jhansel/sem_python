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
    A = self.dof_handler.A[self.k]
    aA1 = self.dof_handler.aA1(U, self.k)
    vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)
    arhoEA = U[self.i_arhoEA]

    arhoABC = vf * self.rho * A
    darhoABC_daA1 = dvf_daA1 * self.rho * A

    arhouABC = vf * self.rho * self.u * A
    darhouABC_daA1 = dvf_daA1 * self.rho * self.u * A

    v, _ = computeSpecificVolume(self.rho)

    E, dE_darhoABC, dE_darhoEA = computeSpecificTotalEnergy(arhoABC, arhoEA)
    dE_daA1 = dE_darhoABC * darhoABC_daA1

    e, _, de_dE = computeSpecificInternalEnergy(self.u, E)
    de_daA1 = de_dE * dE_daA1
    de_darhoEA = de_dE * dE_darhoEA

    p, _, dp_de = self.eos.p(v, e)
    dp_daA1 = dp_de * de_daA1
    dp_darhoEA = dp_de * de_darhoEA

    # momentum
    r[self.i_arhouA] += (arhouABC * self.u + vf * p * A) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhouA,self.i_aA1] += (self.u * darhouABC_daA1 + (dvf_daA1 * p + vf * dp_daA1) * A) * self.nx
    J[self.i_arhouA,self.i_arhoEA] += vf * dp_darhoEA * A * self.nx

    # energy
    r[self.i_arhoEA] += self.u * (arhoEA + vf * p * A) * self.nx
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoEA,self.i_aA1] += self.u * (dvf_daA1 * p + vf * dp_daA1) * A * self.nx
    J[self.i_arhoEA,self.i_arhoEA] += self.u * (1 + vf * dp_darhoEA * A) * self.nx

  def applyStrongBCNonlinearSystem(self, U, r, J):
    A = self.dof_handler.A[self.k]
    aA1 = self.dof_handler.aA1(U, self.k)
    vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)
    arhoA = U[self.i_arhoA]

    arhoABC = vf * self.rho * A
    darhoABC_daA1 = dvf_daA1 * self.rho * A

    # mass
    r[self.i_arhoA] = arhoA - arhoABC
    J[self.i_arhoA,:] = 0
    if (self.model_type == ModelType.TwoPhase):
      J[self.i_arhoA,self.i_aA1] = - darhoABC_daA1
    J[self.i_arhoA,self.i_arhoA] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_arhoA,:] = 0
    A[self.i_arhoA,self.i_arhoA] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    A = self.dof_handler.A[self.k]
    aA1 = self.dof_handler.aA1(U_old, self.k)
    vf, dvf_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)
    arhoA = vf * self.rho * A
    b[self.i_arhoA] = arhoA
