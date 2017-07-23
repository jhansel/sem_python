import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class InterfacePressureParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerParameter("pI_function", "Function for computing interface pressure")

class InterfacePressure(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.pI_function = params.get("pI_function")

  def compute(self, data, der):
    data["pI"], dpI_dp1, dpI_dp2, dpI_dmu = self.pI_function(data["p1"], data["p2"], data["mu"])

    dpI_dvf1 = dpI_dp1 * der["p1"]["vf1"] + dpI_dp2 * der["p2"]["vf1"] + dpI_dmu * der["mu"]["vf1"]
    dpI_darho1 = dpI_dp1 * der["p1"]["arho1"] + dpI_dmu * der["mu"]["arho1"]
    dpI_darho2 = dpI_dp2 * der["p2"]["arho2"] + dpI_dmu * der["mu"]["arho2"]
    dpI_darhou1 = dpI_dp1 * der["p1"]["arhou1"] + dpI_dmu * der["mu"]["arhou1"]
    dpI_darhoE1 = dpI_dp1 * der["p1"]["arhoE1"] + dpI_dmu * der["mu"]["arhoE1"]
    dpI_darhou2 = dpI_dp2 * der["p2"]["arhou2"] + dpI_dmu * der["mu"]["arhou2"]
    dpI_darhoE2 = dpI_dp2 * der["p2"]["arhoE2"] + dpI_dmu * der["mu"]["arhoE2"]
    der["pI"] = {"vf1": dpI_dvf1,
                 "arho1": dpI_darho1, "arhou1": dpI_darhou1, "arhoE1": dpI_darhoE1,
                 "arho2": dpI_darho2, "arhou2": dpI_darhou2, "arhoE2": dpI_darhoE2}
