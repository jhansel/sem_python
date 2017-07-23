import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class ThetaParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerParameter("theta_function", "Function for computing theta")

class Theta(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.theta_function = params.get("theta_function")

  def compute(self, data, der):
    data["theta"], pdtheta_pdvf1, dtheta_dp1, dtheta_dp2 = self.theta_function(
      data["vf1"], data["p1"], data["p2"])

    dtheta_dvf1 = pdtheta_pdvf1 + dtheta_dp1 * der["p1"]["vf1"] + dtheta_dp2 * der["p2"]["vf1"]
    dtheta_darho1 = dtheta_dp1 * der["p1"]["arho1"]
    dtheta_darho2 = dtheta_dp2 * der["p2"]["arho2"]
    dtheta_darhou1 = dtheta_dp1 * der["p1"]["arhou1"]
    dtheta_darhou2 = dtheta_dp2 * der["p2"]["arhou2"]
    dtheta_darhoE1 = dtheta_dp1 * der["p1"]["arhoE1"]
    dtheta_darhoE2 = dtheta_dp2 * der["p2"]["arhoE2"]
    der["theta"] = {"vf1": dtheta_dvf1,
                 "arho1": dtheta_darho1, "arhou1": dtheta_darhou1, "arhoE1": dtheta_darhoE1,
                 "arho2": dtheta_darho2, "arhou2": dtheta_darhou2, "arhoE2": dtheta_darhoE2}
