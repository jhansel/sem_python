from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VolumeFractionPhase2Parameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class VolumeFractionPhase2(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.vf

  def compute(self, data, der):
    data[self.name] = 1.0 - data["vf1"]

    der[self.name] = {"vf1" : 0 * data["vf1"] - 1.0}
