from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from enums import VariableName

class EnergyPressureRelaxationParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class EnergyPressureRelaxation(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoE)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return self.sign * data["pI_bar"] * data["p_relax"] * (data["p1"] - data["p2"]) * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return self.sign * (der["pI_bar"]["vf1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["vf1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * (der["p1"]["vf1"] - der["p2"]["vf1"])) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arho1_index:
      return self.sign * (der["pI_bar"]["arho1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arho1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * der["p1"]["arho1"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhou1_index:
      return self.sign * (der["pI_bar"]["arhou1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhou1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * der["p1"]["arhou1"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoE1_index:
      return self.sign * (der["pI_bar"]["arhoE1"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhoE1"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * data["p_relax"] * der["p1"]["arhoE1"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arho2_index:
      return self.sign * (der["pI_bar"]["arho2"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arho2"] * (data["p1"] - data["p2"]) - \
        data["pI_bar"] * data["p_relax"] * der["p2"]["arho2"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhou2_index:
      return self.sign * (der["pI_bar"]["arhou2"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhou2"] * (data["p1"] - data["p2"]) - \
        data["pI_bar"] * data["p_relax"] * der["p2"]["arhou2"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoE2_index:
      return self.sign * (der["pI_bar"]["arhoE2"] * data["p_relax"] * (data["p1"] - data["p2"]) + \
        data["pI_bar"] * der["p_relax"]["arhoE2"] * (data["p1"] - data["p2"]) - \
        data["pI_bar"] * data["p_relax"] * der["p2"]["arhoE2"]) * \
        data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero