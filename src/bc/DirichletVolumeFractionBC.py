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
  def __init__(self, params, dof_handler, eos_map):
    VolumeFractionBC.__init__(self, params, dof_handler, eos_map)

    self.vf1 = params.get("vf1")

  def apply(self, U, r, J):
    vf1 = U[self.i_vf1]

    r[self.i_vf1] = vf1 - self.vf1
    J[self.i_vf1,:] = 0
    J[self.i_vf1,self.i_vf1] = 1
