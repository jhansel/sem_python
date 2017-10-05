from .OnePhaseBC import OnePhaseBC, OnePhaseBCParameters

class SolidWallBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)

class SolidWallBC(OnePhaseBC):
  def __init__(self, params):
    OnePhaseBC.__init__(self, params)

  def applyWeakBC(self, U, r, J):
    pass

  def applyStrongBCNonlinearSystem(self, U, r, J):
    arhouA = U[self.i_arhouA]

    r[self.i_arhouA] = arhouA
    J[self.i_arhouA,:] = 0
    J[self.i_arhouA,self.i_arhouA] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_arhouA,:] = 0
    A[self.i_arhouA,self.i_arhouA] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    b[self.i_arhouA] = 0
