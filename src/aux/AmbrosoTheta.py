import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoThetaParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerFloatParameter("pressure_relaxation_time", "Relaxation time for pressures")

class AmbrosoTheta(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.pressure_relaxation_time = params.get("pressure_relaxation_time")

  def compute(self, data, der):
    vf1 = data["vf1"]
    p1 = data["p1"]
    p2 = data["p2"]

    vf2 = 1 - vf1
    dvf2_dvf1 = 0 * vf2 - 1

    denominator = self.pressure_relaxation_time * (p1 + p2)
    data["theta"] = vf1 * vf2 / denominator
    dtheta_ddenominator = - vf1 * vf2 / denominator / denominator
    pdtheta_pdvf1 = (vf2 + vf1 * dvf2_dvf1) / denominator
    dtheta_dp1 = dtheta_ddenominator * self.pressure_relaxation_time
    dtheta_dp2 = dtheta_ddenominator * self.pressure_relaxation_time

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
