import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters

sys.path.append(base_dir + "src/base")
from enums import VariableName

class EnergyVolumeFractionGradientParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class EnergyVolumeFractionGradient(Kernel2Phase):
  def __init__(self, params):
    params.set("var_enum", VariableName.ARhoE)
    Kernel2Phase.__init__(self, params)

  def computeResidual(self, data, i):
    return - self.sign * data["pI"] * data["uI"] * data["grad_vf1"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return - self.sign * (der["pI"]["vf1"] * data["uI"] * data["grad_vf1"] * data["phi"][j] \
        + data["pI"] * der["uI"]["vf1"] * data["grad_vf1"] * data["phi"][j] \
        + data["pI"] * data["uI"] * data["grad_phi"][j]) * data["phi"][i] * data["JxW"]
    else:
      aux = - self.sign * data["grad_vf1"] * data["phi"][j] * data["phi"][i] * data["JxW"]
      if var_index == self.arho1_index:
        return (der["pI"]["arho1"] * data["uI"] + data["pI"] * der["uI"]["arho1"]) * aux
      elif var_index == self.arhou1_index:
        return (der["pI"]["arhou1"] * data["uI"] + data["pI"] * der["uI"]["arhou1"]) * aux
      elif var_index == self.arhoE1_index:
        return (der["pI"]["arhoE1"] * data["uI"] + data["pI"] * der["uI"]["arhoE1"]) * aux
      elif var_index == self.arho2_index:
        return (der["pI"]["arho2"] * data["uI"] + data["pI"] * der["uI"]["arho2"]) * aux
      elif var_index == self.arhou2_index:
        return (der["pI"]["arhou2"] * data["uI"] + data["pI"] * der["uI"]["arhou2"]) * aux
      elif var_index == self.arhoE2_index:
        return (der["pI"]["arhoE2"] * data["uI"] + data["pI"] * der["uI"]["arhoE2"]) * aux
      else:
        return self.zero
