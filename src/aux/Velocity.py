import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VelocityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class Velocity(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    arho = data["arho"]
    arhou = data["arhou"]
    data["u"] = arhou / arho

    du_darho = - arhou / arho / arho
    du_darhou = 1.0 / arho
    der["u"] = {"arho" : du_darho, "arhou" : du_darhou}
