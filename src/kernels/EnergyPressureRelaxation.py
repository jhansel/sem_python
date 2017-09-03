from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from enums import VariableName

class EnergyPressureRelaxationParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class EnergyPressureRelaxation(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoEA)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return self.sign * data["pI_bar"] * data["p_relax"] * (data["p1"] - data["p2"]) * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return self.sign * (der["pI_bar"]["vf1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["vf1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * (der["p1"]["vf1"] - der["p2"]["vf1"])) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoA1_index:
      return self.sign * (der["pI_bar"]["arhoA1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhoA1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * der["p1"]["arhoA1"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhouA1_index:
      return self.sign * (der["pI_bar"]["arhouA1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhouA1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * der["p1"]["arhouA1"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoEA1_index:
      return self.sign * (der["pI_bar"]["arhoEA1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhoEA1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * der["p1"]["arhoEA1"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoA2_index:
      return self.sign * (der["pI_bar"]["arhoA2"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhoA2"] * (data["p1"] - data["p2"]) - \
        data["pI_bar"] * data["p_relax"] * der["p2"]["arhoA2"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhouA2_index:
      return self.sign * (der["pI_bar"]["arhouA2"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhouA2"] * (data["p1"] - data["p2"]) - \
        data["pI_bar"] * data["p_relax"] * der["p2"]["arhouA2"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoEA2_index:
      return self.sign * (der["pI_bar"]["arhoEA2"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhoEA2"] * (data["p1"] - data["p2"]) - \
        data["pI_bar"] * data["p_relax"] * der["p2"]["arhoEA2"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
