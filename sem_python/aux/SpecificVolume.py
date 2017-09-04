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
    dv_dvf1 = dv_drho * der[self.rho]["vf1"]
    dv_darho = dv_drho * der[self.rho][self.arho]
    der[self.name]["vf1"] = dv_dvf1
    der[self.name][self.arho] = dv_darho
