import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class InternalEnergyDensityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

## Density multiplied by specific internal energy: \f$\rho e\f$
class InternalEnergyDensity(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    data[self.rhoe] = data[self.rho] * data[self.e]

    drhoe_dvf1 = der[self.rho]["vf1"] * data[self.e]
    drhoe_darho = der[self.rho][self.arho] * data[self.e] + data[self.rho] * der[self.e][self.arho]
    drhoe_darhou = data[self.rho] * der[self.e][self.arhou]
    drhoe_darhoE = data[self.rho] * der[self.e][self.arhoE]

    der[self.rhoe] = {"vf1": drhoe_dvf1, self.arho: drhoe_darho, self.arhou: drhoe_darhou, self.arhoE: drhoe_darhoE}
