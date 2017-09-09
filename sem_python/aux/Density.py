from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class DensityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class Density(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.rho

  def compute(self, data, der):
    vf = data[self.vf]
    data[self.name] = data[self.arhoA] / vf

    drho_dvf = - data[self.arhoA] / vf / vf
    drho_daA1 = drho_dvf * der[self.vf]["aA1"]
    drho_darhoA = 1.0 / vf
    der[self.name]["aA1"] = drho_daA1
    der[self.name][self.arhoA] = drho_darhoA
