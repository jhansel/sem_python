from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VolumeFractionPhase1Parameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class VolumeFractionPhase1(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.vf

  def compute(self, data, der):
    data[self.name] = data["aA1"] / data["A"]

    # The derivative w.r.t. area is needed because this derivative will be
    # used by AuxGradient to compute the gradient of volume fraction.
    der[self.name]["A"] = - data["aA1"] / data["A"]**2
    der[self.name]["aA1"] = 1.0 / data["A"]
