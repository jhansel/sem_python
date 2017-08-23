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
    dp_darho = dp_dv * der[self.v][self.arho] + dp_de * der[self.e][self.arho]
    dp_darhou = dp_de * der[self.e][self.arhou]
    dp_darhoE = dp_de * der[self.e][self.arhoE]
    der[self.name] = {"vf1" : dp_dvf1, self.arho : dp_darho, self.arhou : dp_darhou, self.arhoE : dp_darhoE}
