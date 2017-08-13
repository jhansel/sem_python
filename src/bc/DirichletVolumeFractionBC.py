import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/bc")
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
    vf1 = U[self.i_vf1]

    r[self.i_vf1] = vf1 - self.vf1
    J[self.i_vf1,:] = 0
    J[self.i_vf1,self.i_vf1] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_vf1,:] = 0
    A[self.i_vf1,self.i_vf1] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    b[self.i_vf1] = self.vf1
