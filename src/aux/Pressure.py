from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class PressureParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)
    self.registerFunctionParameter("p_function", "Pressure function from EOS")

class Pressure(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.p
    self.p_function = params.get("p_function")

  def compute(self, data, der):
    p, dp_dv, dp_de = self.p_function(data[self.v], data[self.e])
    data[self.name] = p

    dp_dvf1 = dp_dv * der[self.v]["vf1"]
    dp_darhoA = dp_dv * der[self.v][self.arhoA] + dp_de * der[self.e][self.arhoA]
    dp_darhouA = dp_de * der[self.e][self.arhouA]
    dp_darhoEA = dp_de * der[self.e][self.arhoEA]
    der[self.name]["vf1"] = dp_dvf1
    der[self.name][self.arhoA] = dp_darhoA
    der[self.name][self.arhouA] = dp_darhouA
    der[self.name][self.arhoEA] = dp_darhoEA
