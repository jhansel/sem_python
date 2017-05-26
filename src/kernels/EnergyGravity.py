import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters

sys.path.append(base_dir + "src/base")
from enums import VariableName

class EnergyGravityParameters(Kernel1PhaseParameters):
  def __init__(self):
    Kernel1PhaseParameters.__init__(self)

class EnergyGravity(Kernel1Phase):
  def __init__(self, params, dof_handler):
    Kernel1Phase.__init__(self, params, dof_handler, VariableName.ARhoE)

  def computeResidual(self, data, i):
    return -data[self.arhou] * data["g"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.arhou_index:
      return -data["g"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
