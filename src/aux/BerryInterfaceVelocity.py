from numpy import sign

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfaceVelocityParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryInterfaceVelocity(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    u_int_bar = data["u_int_bar"]
    grad_vf1 = data["grad_vf1"]
    p1 = data["p1"]
    p2 = data["p2"]
    z1 = data["z1"]
    z2 = data["z2"]

    numerator = sign(grad_vf1) * (p2 - p1)
    denominator = z1 + z2
    data["u_int"] = u_int_bar + numerator / denominator

    dnumerator = 1.0 / denominator
    ddenominator = -numerator / denominator**2
    dvf1 = der["u_int_bar"]["vf1"] + dnumerator * sign(grad_vf1) * (der["p2"]["vf1"] - der["p1"]["vf1"]) \
      + ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    darho1 = der["u_int_bar"]["arho1"] - dnumerator * sign(grad_vf1) * der["p1"]["arho1"] + ddenominator * der["z1"]["arho1"]
    darhou1 = der["u_int_bar"]["arhou1"] - dnumerator * sign(grad_vf1) * der["p1"]["arhou1"] + ddenominator * der["z1"]["arhou1"]
    darhoE1 = der["u_int_bar"]["arhoE1"] - dnumerator * sign(grad_vf1) * der["p1"]["arhoE1"] + ddenominator * der["z1"]["arhoE1"]
    darho2 = der["u_int_bar"]["arho2"] + dnumerator * sign(grad_vf1) * der["p2"]["arho2"] + ddenominator * der["z2"]["arho2"]
    darhou2 = der["u_int_bar"]["arhou2"] + dnumerator * sign(grad_vf1) * der["p2"]["arhou2"] + ddenominator * der["z2"]["arhou2"]
    darhoE2 = der["u_int_bar"]["arhoE2"] + dnumerator * sign(grad_vf1) * der["p2"]["arhoE2"] + ddenominator * der["z2"]["arhoE2"]

    der["u_int"] = {"vf1": dvf1, "arho1": darho1, "arhou1": darhou1, "arhoE1": darhoE1,
      "arho2": darho2, "arhou2": darhou2, "arhoE2": darhoE2}
