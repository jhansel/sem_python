from Kernel1Phase import Kernel1Phase, Kernel1PhaseParameters
from enums import VariableName

class EnergyHeatTransferParameters(Kernel1PhaseParameters):
  def __init__(self):
    Kernel1PhaseParameters.__init__(self)

class EnergyHeatTransfer(Kernel1Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoEA)
    Kernel1Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return -data["htc_wall"] * (data["T_wall"] - data[self.T]) * data["P_heat"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.aA1_index:
      return data["htc_wall"] * der[self.T]["aA1"] * data["P_heat"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoA_index:
      return data["htc_wall"] * der[self.T][self.arhoA] * data["P_heat"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhouA_index:
      return data["htc_wall"] * der[self.T][self.arhouA] * data["P_heat"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoEA_index:
      return data["htc_wall"] * der[self.T][self.arhoEA] * data["P_heat"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
