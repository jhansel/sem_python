from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from enums import VariableName

class EnergyVolumeFractionGradientParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class EnergyVolumeFractionGradient(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoEA)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return - self.sign * data["pI"] * data["uI"] * data["grad_vf1"] * data["A"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.aA1_index:
      return - self.sign * (der["pI"]["aA1"] * data["uI"] * data["grad_vf1"] * data["phi"][j] \
        + data["pI"] * der["uI"]["aA1"] * data["grad_vf1"] * data["phi"][j] \
        + data["pI"] * data["uI"] * der["vf1"]["aA1"] * data["grad_phi"][j]) * data["A"] * data["phi"][i] * data["JxW"]
    else:
      aux = - self.sign * data["grad_vf1"] * data["A"] * data["phi"][j] * data["phi"][i] * data["JxW"]
      if var_index == self.arhoA1_index:
        return (der["pI"]["arhoA1"] * data["uI"] + data["pI"] * der["uI"]["arhoA1"]) * aux
      elif var_index == self.arhouA1_index:
        return (der["pI"]["arhouA1"] * data["uI"] + data["pI"] * der["uI"]["arhouA1"]) * aux
      elif var_index == self.arhoEA1_index:
        return (der["pI"]["arhoEA1"] * data["uI"] + data["pI"] * der["uI"]["arhoEA1"]) * aux
      elif var_index == self.arhoA2_index:
        return (der["pI"]["arhoA2"] * data["uI"] + data["pI"] * der["uI"]["arhoA2"]) * aux
      elif var_index == self.arhouA2_index:
        return (der["pI"]["arhouA2"] * data["uI"] + data["pI"] * der["uI"]["arhouA2"]) * aux
      elif var_index == self.arhoEA2_index:
        return (der["pI"]["arhoEA2"] * data["uI"] + data["pI"] * der["uI"]["arhoEA2"]) * aux
      else:
        return self.zero
