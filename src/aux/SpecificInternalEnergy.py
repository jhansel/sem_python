import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SpecificInternalEnergyParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class SpecificInternalEnergy(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    u = data["u"]
    E = data["E"]
    data["e"] = E - 0.5 * u * u

    de_dE = 1.0
    de_du = - u
    de_darho = de_dE * der["E"]["arho"] + de_du * der["u"]["arho"]
    de_darhou = de_du * der["u"]["arhou"]
    de_darhoE = de_dE * der["E"]["arhoE"]
    der["e"] = {"arho" : de_darho, "arhou" : de_darhou, "arhoE" : de_darhoE}
