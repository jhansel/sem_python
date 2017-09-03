from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class InternalEnergyDensityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

## Density multiplied by specific internal energy: \f$\rho e\f$
class InternalEnergyDensity(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.rhoe

  def compute(self, data, der):
    data[self.name] = data[self.rho] * data[self.e]

    drhoe_dvf1 = der[self.rho]["vf1"] * data[self.e]
    drhoe_darho = der[self.rho][self.arho] * data[self.e] + data[self.rho] * der[self.e][self.arho]
    drhoe_darhou = data[self.rho] * der[self.e][self.arhou]
    drhoe_darhoE = data[self.rho] * der[self.e][self.arhoE]

    der[self.name]["vf1"] = drhoe_dvf1
    der[self.name][self.arho] = drhoe_darho
    der[self.name][self.arhou] = drhoe_darhou
    der[self.name][self.arhoE] = drhoe_darhoE
