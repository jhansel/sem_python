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
    data[self.name] = data[self.arho] / vf

    drho_dvf = - data[self.arho] / vf / vf
    drho_dvf1 = drho_dvf * der[self.vf]["vf1"]
    drho_darho = 1.0 / vf
    der[self.name] = {"vf1" : drho_dvf1, self.arho : drho_darho}
