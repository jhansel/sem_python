from VolumeFractionBC import VolumeFractionBC, VolumeFractionBCParameters

class DirichletVolumeFractionBCParameters(VolumeFractionBCParameters):
  def __init__(self):
    VolumeFractionBCParameters.__init__(self)
    self.registerFloatParameter("vf1", "Specified volume fraction of first phase")

class DirichletVolumeFractionBC(VolumeFractionBC):
  def __init__(self, params):
    VolumeFractionBC.__init__(self, params)

    self.vf1 = params.get("vf1")

  def applyWeakBC(self, U, r, J):
    pass

  def applyStrongBCNonlinearSystem(self, U, r, J):
    A = self.dof_handler.A[self.k]
    aA1 = self.dof_handler.aA1(U, self.k)

    vf1, dvf1_daA1 = computeVolumeFraction(aA1, A, self.phase, self.model_type)

    r[self.i_aA1] = vf1 - self.vf1
    J[self.i_aA1,:] = 0
    J[self.i_aA1,self.i_aA1] = dvf1_daA1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_aA1,:] = 0
    A[self.i_aA1,self.i_aA1] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    A = self.dof_handler.A[self.k]

    b[self.i_aA1] = self.vf1 * A
