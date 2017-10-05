from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SpecificVolumeParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class SpecificVolume(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.v

  def compute(self, data, der):
    rho = data[self.rho]
    data[self.name] = 1.0 / rho

    dv_drho = - 1.0 / rho / rho
    dv_daA1 = dv_drho * der[self.rho]["aA1"]
    dv_darhoA = dv_drho * der[self.rho][self.arhoA]
    der[self.name]["aA1"] = dv_daA1
    der[self.name][self.arhoA] = dv_darhoA
