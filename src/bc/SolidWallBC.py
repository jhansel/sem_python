import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/bc")
from OnePhaseBC import OnePhaseBC, OnePhaseBCParameters

class SolidWallBCParameters(OnePhaseBCParameters):
  def __init__(self):
    OnePhaseBCParameters.__init__(self)

class SolidWallBC(OnePhaseBC):
  def __init__(self, params, dof_handler, eos_map):
    OnePhaseBC.__init__(self, params, dof_handler, eos_map)

  def applyWeakBC(self, U, r, J):
    pass

  def applyStrongBCNonlinearSystem(self, U, r, J):
    arhou = U[self.i_arhou]

    r[self.i_arhou] = arhou
    J[self.i_arhou,:] = 0
    J[self.i_arhou,self.i_arhou] = 1

  def applyStrongBCLinearSystemMatrix(self, A):
    A[self.i_arhou,:] = 0
    A[self.i_arhou,self.i_arhou] = 1

  def applyStrongBCLinearSystemRHSVector(self, U_old, b):
    b[self.i_arhou] = 0
