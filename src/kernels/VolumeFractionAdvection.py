import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from Kernel2Phase import Kernel2Phase, Kernel2PhaseParameters

sys.path.append(base_dir + "src/base")
from enums import VariableName

class VolumeFractionAdvectionParameters(Kernel2PhaseParameters):
  def __init__(self):
    Kernel2PhaseParameters.__init__(self)

class VolumeFractionAdvection(Kernel2Phase):
  def __init__(self, params, dof_handler):
    Kernel2Phase.__init__(self, params, dof_handler, VariableName.VF1)

  def computeResidual(self, data, i):
    return (data["uI"] * data["dvf1_dx"] - data["theta"] * (data["p1"] - data["p2"])) \
      * data["phi"][i] * data["JxW"]

  def computeJacobian(self, data, der, var_index, i, j):
    if var_index == self.vf1_index:
      return (data["uI"] * data["grad_phi"][j] \
        - der["theta"]["vf1"] * (data["p1"] - data["p2"]) * data["phi"][j] \
        - data["theta"] * (der["p1"]["vf1"] - der["p2"]["vf1"]) * data["phi"][j]) * data["phi"][i] * data["JxW"]
    elif var_index == self.arho1_index:
      return (der["uI"]["arho1"] * data["dvf1_dx"] - der["theta"]["arho1"] * (data["p1"] - data["p2"]) \
        - data["theta"] * der["p1"]["arho1"]) * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhou1_index:
      return (der["uI"]["arhou1"] * data["dvf1_dx"] - der["theta"]["arhou1"] * (data["p1"] - data["p2"]) \
        - data["theta"] * der["p1"]["arhou1"]) * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoE1_index:
      return (- der["theta"]["arhoE1"] * (data["p1"] - data["p2"]) \
        - data["theta"] * der["p1"]["arhoE1"]) * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arho2_index:
      return (der["uI"]["arho2"] * data["dvf1_dx"] - der["theta"]["arho2"] * (data["p1"] - data["p2"]) \
        + data["theta"] * der["p2"]["arho2"]) * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhou2_index:
      return (der["uI"]["arhou2"] * data["dvf1_dx"] - der["theta"]["arhou2"] * (data["p1"] - data["p2"]) \
        + data["theta"] * der["p2"]["arhou2"]) * data["phi"][j] * data["phi"][i] * data["JxW"]
    elif var_index == self.arhoE2_index:
      return (- der["theta"]["arhoE2"] * (data["p1"] - data["p2"]) \
        + data["theta"] * der["p2"]["arhoE2"]) * data["phi"][j] * data["phi"][i] * data["JxW"]
    else:
      return self.zero
