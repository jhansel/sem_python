from Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters
from enums import VariableName

class MomentumAdvectionParameters(Kernel1PhaseParameters):
  def __init__(self):
    Kernel1PhaseParameters.__init__(self)

class MomentumAdvection(Kernel1Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoU)
    Kernel1Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return -(data[self.arhou] * data[self.u] + data[self.vf] * data[self.p]) * data["grad_phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return -(der[self.vf]["vf1"] * data[self.p] \
        + data[self.vf] * der[self.p]["vf1"]) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arho_index:
      return -(data[self.arhou] * der[self.u][self.arho] \
        + data[self.vf] * der[self.p][self.arho]) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhou_index:
      return -(data[self.u] + data[self.arhou] * der[self.u][self.arhou] \
        + data[self.vf] * der[self.p][self.arhou]) * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    elif var_index == self.arhoE_index:
      return -data[self.vf] * der[self.p][self.arhoE] * data["phi"][j] * data["grad_phi"][i] * data["JxW"]
    else:
      return self.zero
