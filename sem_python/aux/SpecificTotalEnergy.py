from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SpecificTotalEnergyParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class SpecificTotalEnergy(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.E

  def compute(self, data, der):
    arho = data[self.arho]
    arhoE = data[self.arhoE]
    data[self.name] = arhoE / arho

    dE_darho = - arhoE / arho / arho
    dE_darhoE = 1.0 / arho
    der[self.name][self.arho] = dE_darho
    der[self.name][self.arhoE] = dE_darhoE
