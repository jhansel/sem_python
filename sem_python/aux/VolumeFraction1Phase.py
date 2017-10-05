from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VolumeFraction1PhaseParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class VolumeFraction1Phase(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.vf

  def compute(self, data, der):
    data[self.name] = 1

    der[self.name]["A"] = 0
    der[self.name]["aA1"] = 0
