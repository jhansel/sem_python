from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class SoundSpeedParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)
    self.registerFunctionParameter("c_function", "Sound speed function from EOS")

class SoundSpeed(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.c
    self.c_function = params.get("c_function")

  def compute(self, data, der):
    c, dc_dv, dc_dp = self.c_function(data[self.v], data[self.p])
    data[self.name] = c

    dc_daA1 = dc_dv * der[self.v]["aA1"] + dc_dp * der[self.p]["aA1"]
    dc_darhoA = dc_dv * der[self.v][self.arhoA] + dc_dp * der[self.p][self.arhoA]
    dc_darhouA = dc_dp * der[self.p][self.arhouA]
    dc_darhoEA = dc_dp * der[self.p][self.arhoEA]
    der[self.name]["aA1"] = dc_daA1
    der[self.name][self.arhoA] = dc_darhoA
    der[self.name][self.arhouA] = dc_darhouA
    der[self.name][self.arhoEA] = dc_darhoEA
