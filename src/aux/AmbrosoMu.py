import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoMuParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoMu(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    beta = data["beta"]
    T1 = data["T1"]
    T2 = data["T2"]

    denominator = beta * T1 + (1 - beta) * T2
    data["mu"] = (1 - beta) * T2 / denominator
    dmu_dT1 = - (1 - beta) * T2 / denominator / denominator * beta
    dmu_dT2 = (1 - beta) / denominator - (1 - beta) * T2 / denominator / denominator * (1 - beta)
    dmu_dbeta = - T2 / denominator - (1 - beta) * T2 / denominator / denominator * (T1 - T2)

    dmu_dvf1 = dmu_dT1 * der["T1"]["vf1"] + dmu_dT2 * der["T2"]["vf1"]
    dmu_darho1 = dmu_dT1 * der["T1"]["arho1"] + dmu_dbeta * der["beta"]["arho1"]
    dmu_darho2 = dmu_dT2 * der["T2"]["arho2"] + dmu_dbeta * der["beta"]["arho2"]
    dmu_darhou1 = dmu_dT1 * der["T1"]["arhou1"]
    dmu_darhou2 = dmu_dT2 * der["T2"]["arhou2"]
    dmu_darhoE1 = dmu_dT1 * der["T1"]["arhoE1"]
    dmu_darhoE2 = dmu_dT2 * der["T2"]["arhoE2"]

    der["mu"] = {"vf1": dmu_dvf1,
                 "arho1": dmu_darho1, "arhou1": dmu_darhou1, "arhoE1": dmu_darhoE1,
                 "arho2": dmu_darho2, "arhou2": dmu_darhou2, "arhoE2": dmu_darhoE2}
