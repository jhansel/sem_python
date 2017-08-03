import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfaceVelocityBarParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryInterfaceVelocityBar(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    u1 = data["u1"]
    u2 = data["u2"]
    z1 = data["z1"]
    z2 = data["z2"]

    numerator = (z1 * u1 + z2 * u2)
    data["uI_bar"] = numerator / (z1 + z2)

    ddenominator = -1.0 / (z1 + z2)**2
    dvf1 = (der["z1"]["vf1"] * u1 + der["z2"]["vf1"] * u2) / (z1 + z2) \
      + numerator * ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    darho1 = (der["z1"]["arho1"] * u1 + z1 * der["u1"]["arho1"]) / (z1 + z2) \
      + numerator * ddenominator * der["z1"]["arho1"]
    darhou1 = (der["z1"]["arhou1"] * u1 + z1 * der["u1"]["arhou1"]) / (z1 + z2) \
      + numerator * ddenominator * der["z1"]["arhou1"]
    darhoE1 = der["z1"]["arhoE1"] * u1 / (z1 + z2) + numerator * ddenominator * der["z1"]["arhoE1"]
    darho2 = (der["z2"]["arho2"] * u2 + z2 * der["u2"]["arho2"]) / (z1 + z2) \
      + numerator * ddenominator * der["z2"]["arho2"]
    darhou2 = (der["z2"]["arhou2"] * u2 + z2 * der["u2"]["arhou2"]) / (z1 + z2) \
      + numerator * ddenominator * der["z2"]["arhou2"]
    darhoE2 = der["z2"]["arhoE2"] * u2 / (z1 + z2) + numerator * ddenominator * der["z2"]["arhoE2"]

    der["uI_bar"] = {"vf1": dvf1, "arho1": darho1, "arhou1": darhou1, "arhoE1": darhoE1,
      "arho2": darho2, "arhou2": darhou2, "arhoE2": darhoE2}
