from .Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters
from ..base.enums import VariableName

class MomentumGravityParameters(Kernel1PhaseParameters):
  def __init__(self):
    Kernel1PhaseParameters.__init__(self)

class MomentumGravity(Kernel1Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoU)
    Kernel1Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return -data[self.arho] * data["g"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.arho_index:
      return -data["g"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
