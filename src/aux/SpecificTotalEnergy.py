from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SpecificTotalEnergyParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class SpecificTotalEnergy(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    arho = data[self.arho]
    arhoE = data[self.arhoE]
    data[self.E] = arhoE / arho

    dE_darho = - arhoE / arho / arho
    dE_darhoE = 1.0 / arho
    der[self.E] = {self.arho : dE_darho, self.arhoE : dE_darhoE}
