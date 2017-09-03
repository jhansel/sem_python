from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from enums import VariableName

class MomentumVolumeFractionGradientParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class MomentumVolumeFractionGradient(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoUA)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return - self.sign * data["pI"] * data["grad_vf1"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return (-self.sign * der["pI"]["vf1"] * data["grad_vf1"] * data["phi"][j] \
        - self.sign * data["pI"] * data["grad_phi"][j]) * data["phi"][i] * data["JxW"]
    else:
      aux = - self.sign * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
      if var_index == self.arhoA1_index:
        return der["pI"]["arhoA1"]  * aux
      elif var_index == self.arhouA1_index:
        return der["pI"]["arhouA1"] * aux
      elif var_index == self.arhoEA1_index:
        return der["pI"]["arhoEA1"] * aux
      elif var_index == self.arhoA2_index:
        return der["pI"]["arhoA2"]  * aux
      elif var_index == self.arhouA2_index:
        return der["pI"]["arhouA2"] * aux
      elif var_index == self.arhoEA2_index:
        return der["pI"]["arhoEA2"] * aux
      else:
        return self.zero
