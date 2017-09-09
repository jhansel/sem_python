from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from enums import VariableName

class VolumeFractionPressureRelaxationParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class VolumeFractionPressureRelaxation(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.AA1)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return -data["p_relax"] * (data["p1"] - data["p2"]) * data["A"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.aA1_index:
      return -(data["p_relax"] * (der["p1"]["aA1"] - der["p2"]["aA1"]) + \
        der["p_relax"]["aA1"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoA1_index:
      return -(data["p_relax"] * der["p1"]["arhoA1"] + \
        der["p_relax"]["arhoA1"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhouA1_index:
      return -(data["p_relax"] * der["p1"]["arhouA1"] + \
        der["p_relax"]["arhouA1"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoEA1_index:\
      return -(data["p_relax"] * der["p1"]["arhoEA1"] + \
        der["p_relax"]["arhoEA1"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoA2_index:
      return -(-data["p_relax"] * der["p2"]["arhoA2"] + \
        der["p_relax"]["arhoA2"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhouA2_index:
      return -(-data["p_relax"] * der["p2"]["arhouA2"] + \
        der["p_relax"]["arhouA2"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoEA2_index:
      return -(-data["p_relax"] * der["p2"]["arhoEA2"] + \
        der["p_relax"]["arhoEA2"] * (data["p1"] - data["p2"])) * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
