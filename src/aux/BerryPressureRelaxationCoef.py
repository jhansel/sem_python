import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryPressureRelaxationCoefParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryPressureRelaxationCoef(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    a_int = data["a_int"]
    z1 = data["z1"]
    z2 = data["z2"]

    data["p_relax"] = a_int / (z1 + z2)

    ddenominator = - a_int / (z1 + z2)**2
    dvf1 = der["a_int"]["vf1"] / (z1 + z2) + ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    darho1 = ddenominator * der["z1"]["arho1"]
    darhou1 = ddenominator * der["z1"]["arhou1"]
    darhoE1 = ddenominator * der["z1"]["arhoE1"]
    darho2 = ddenominator * der["z2"]["arho2"]
    darhou2 = ddenominator * der["z2"]["arhou2"]
    darhoE2 = ddenominator * der["z2"]["arhoE2"]

    der["p_relax"] = {"vf1": dvf1, "arho1": darho1, "arhou1": darhou1, "arhoE1": darhoE1,
      "arho2": darho2, "arhou2": darhou2, "arhoE2": darhoE2}
