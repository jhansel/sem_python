import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class DensityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class Density(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    vf = data["vf"]
    data["rho"] = data["arho"] / vf

    drho_dvf = - data["arho"] / vf / vf
    drho_dvf1 = drho_dvf * der["vf"]["vf1"]
    drho_darho = 1.0 / vf
    der["rho"] = {"vf1" : drho_dvf1, "arho" : drho_darho}
