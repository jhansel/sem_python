import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import VariableName

sys.path.append(base_dir + "src/bc")
from BC import BC, BCParameters

class VolumeFractionBCParameters(BCParameters):
  def __init__(self):
    BCParameters.__init__(self)

class VolumeFractionBC(BC):
  def __init__(self, params, dof_handler, eos_map):
    BC.__init__(self, params, dof_handler, eos_map)

    self.i_vf1 = self.dof_handler.i(self.k, VariableName.VF1)
