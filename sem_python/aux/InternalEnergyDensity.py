from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

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

    drhoe_daA1 = der[self.rho]["aA1"] * data[self.e]
    drhoe_darhoA = der[self.rho][self.arhoA] * data[self.e] + data[self.rho] * der[self.e][self.arhoA]
    drhoe_darhouA = data[self.rho] * der[self.e][self.arhouA]
    drhoe_darhoEA = data[self.rho] * der[self.e][self.arhoEA]

    der[self.name]["aA1"] = drhoe_daA1
    der[self.name][self.arhoA] = drhoe_darhoA
    der[self.name][self.arhouA] = drhoe_darhouA
    der[self.name][self.arhoEA] = drhoe_darhoEA
