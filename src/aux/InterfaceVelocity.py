import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class InterfaceVelocityParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerParameter("uI_function", "Function for computing interface velocity")

class InterfaceVelocity(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.uI_function = params.get("uI_function")

  def compute(self, data, der):
    data["uI"], duI_du1, duI_du2, duI_dbeta = self.uI_function(data["u1"], data["u2"], data["beta"])

    duI_darho1 = duI_du1 * der["u1"]["arho1"] + duI_dbeta * der["beta"]["arho1"]
    duI_darho2 = duI_du2 * der["u2"]["arho2"] + duI_dbeta * der["beta"]["arho2"]
    duI_darhou1 = duI_du1 * der["u1"]["arhou1"]
    duI_darhou2 = duI_du2 * der["u2"]["arhou2"]
    der["uI"] = {"arho1": duI_darho1, "arho2": duI_darho2, "arhou1": duI_darhou1, "arhou2": duI_darhou2}
