from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SpecificInternalEnergyParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class SpecificInternalEnergy(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    u = data[self.u]
    E = data[self.E]
    data[self.e] = E - 0.5 * u * u

    de_dE = 1.0
    de_du = - u
    de_darho = de_dE * der[self.E][self.arho] + de_du * der[self.u][self.arho]
    de_darhou = de_du * der[self.u][self.arhou]
    de_darhoE = de_dE * der[self.E][self.arhoE]
    der[self.e] = {self.arho : de_darho, self.arhou : de_darhou, self.arhoE : de_darhoE}
