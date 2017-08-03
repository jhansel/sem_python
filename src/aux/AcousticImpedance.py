import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class AcousticImpedanceParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class AcousticImpedance(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    data[self.rhoc] = data[self.rho] * data[self.c]

    drhoc_dvf1 = der[self.rho]["vf1"] * data[self.c] + data[self.rho] * der[self.c]["vf1"]
    drhoc_darho = der[self.rho][self.arho] * data[self.c] + data[self.rho] * der[self.c][self.arho]
    drhoc_darhou = data[self.rho] * der[self.c][self.arhou]
    drhoc_darhoE = data[self.rho] * der[self.c][self.arhoE]

    der[self.rhoc] = {"vf1": drhoc_dvf1, self.arho: drhoc_darho, self.arhou: drhoc_darhou, self.arhoE: drhoc_darhoE}
