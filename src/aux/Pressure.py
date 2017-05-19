import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class PressureParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class Pressure(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    p, dp_dv, dp_de = self.p_function(data["v"], data["e"])
    data["p"] = p

    dp_dvf1 = dp_dv * der["v"]["vf1"]
    dp_darho = dp_dv * der["v"]["arho"] + dp_de * der["e"]["arho"]
    dp_darhou = dp_de * der["e"]["arhou"]
    dp_darhoE = dp_de * der["e"]["arhoE"]
    der["p"] = {"vf1" : dp_dvf1, "arho" : dp_darho, "arhou" : dp_darhou, "arhoE" : dp_darhoE}
