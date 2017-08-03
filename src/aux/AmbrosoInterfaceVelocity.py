import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoInterfaceVelocityParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoInterfaceVelocity(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    beta = data["beta"]
    u1 = data["u1"]
    u2 = data["u2"]

    data["uI"] = beta * u1 + (1 - beta) * u2
    duI_du1 = beta
    duI_du2 = (1 - beta)
    duI_dbeta = u1 - u2

    duI_darho1 = duI_du1 * der["u1"]["arho1"] + duI_dbeta * der["beta"]["arho1"]
    duI_darho2 = duI_du2 * der["u2"]["arho2"] + duI_dbeta * der["beta"]["arho2"]
    duI_darhou1 = duI_du1 * der["u1"]["arhou1"]
    duI_darhou2 = duI_du2 * der["u2"]["arhou2"]
    der["uI"] = {"arho1": duI_darho1, "arho2": duI_darho2, "arhou1": duI_darhou1, "arhou2": duI_darhou2}
