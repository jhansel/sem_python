import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters

sys.path.append(base_dir + "src/base")
from enums import VariableName

class MomentumVolumeFractionGradientParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class MomentumVolumeFractionGradient(Kernel2Phase):
  def __init__(self, params, dof_handler):
    Kernel2Phase.__init__(self, params, dof_handler, VariableName.ARhoU)

  def computeResidual(self, data, i):
    return - self.sign * data["pI"] * data["dvf1_dx"] * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return (-self.sign * der["pI"]["vf1"] * data["dvf1_dx"] * data["phi"][j] \
        - self.sign * data["pI"] * data["grad_phi"][j]) * data["phi"][i] * data["JxW"]
    else:
      aux = - self.sign * data["dvf1_dx"] * data["phi"][j] * data["phi"][i] * data["JxW"]
      if var_index == self.arho1_index:
        return der["pI"]["arho1"]  * aux
      elif var_index == self.arhou1_index:
        return der["pI"]["arhou1"] * aux
      elif var_index == self.arhoE1_index:
        return der["pI"]["arhoE1"] * aux
      elif var_index == self.arho2_index:
        return der["pI"]["arho2"]  * aux
      elif var_index == self.arhou2_index:
        return der["pI"]["arhou2"] * aux
      elif var_index == self.arhoE2_index:
        return der["pI"]["arhoE2"] * aux
      else:
        return self.zero
