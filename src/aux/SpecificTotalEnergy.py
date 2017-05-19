import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SpecificTotalEnergyParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class SpecificTotalEnergy(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    arho = data["arho"]
    arhoE = data["arhoE"]
    data["E"] = arhoE / arho

    dE_darho = - arhoE / arho / arho
    dE_darhoE = 1.0 / arho
    der["E"] = {"arho" : dE_darho, "arhoE" : dE_darhoE}
