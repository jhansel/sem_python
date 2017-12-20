from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class MassFluxParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class MassFlux(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.inviscflux_arhoA

  def compute(self, data, der):
    data[self.name] = data[self.arhouA]
    der[self.name][self.arhouA] = 1.0
