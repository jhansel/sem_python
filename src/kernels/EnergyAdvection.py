import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters

sys.path.append(base_dir + "src/base")
from enums import VariableName

class EnergyAdvectionParameters(Kernel1PhaseParameters):
  def __init__(self):
    Kernel1PhaseParameters.__init__(self)

class EnergyAdvection(Kernel1Phase):
  def __init__(self, params, dof_handler):
    Kernel1Phase.__init__(self, params, dof_handler, VariableName.ARhoE)

  def computeResidual(self, data, i):
    return -data[self.u] * (data[self.arhoE] + data[self.vf] * data[self.p]) * data["grad_phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return -data[self.u] * (data[self.vf] * der[self.p]["vf1"] \
        + der[self.vf]["vf1"] * data[self.p]) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arho_index:
      return -(data[self.u] * (data[self.vf] * der[self.p][self.arho]) \
        + der[self.u][self.arho] * (data[self.arhoE] + data[self.vf] * data[self.p])) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhou_index:
      return -(data[self.u] * (data[self.vf] * der[self.p][self.arhou]) \
        + der[self.u][self.arhou] * (data[self.arhoE] + data[self.vf] * data[self.p])) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhoE_index:
      return -data[self.u] * (1 + data[self.vf] * der[self.p][self.arhoE]) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    else:
      return self.zero
