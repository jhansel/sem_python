from .Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters
from ..base.enums import VariableName

class VolumeFractionAdvectionParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class VolumeFractionAdvection(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.VF1)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return data["uI"] * data["grad_vf1"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return (der["uI"]["vf1"] * data["grad_vf1"] * data["phi"][j] + data["uI"] * data["grad_phi"][j]) * data["phi"][i] * data["JxW"]
    elif var_index == self.arho1_index:
      return der["uI"]["arho1"] * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhou1_index:
      return der["uI"]["arhou1"] * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoE1_index:
      return der["uI"]["arhoE1"] * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arho2_index:
      return der["uI"]["arho2"] * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhou2_index:
      return der["uI"]["arhou2"] * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoE2_index:
      return der["uI"]["arhoE2"] * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
