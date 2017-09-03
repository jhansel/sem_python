from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryPressureRelaxationCoefParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryPressureRelaxationCoef(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "p_relax"

  def compute(self, data, der):
    a_int = data["a_int"]
    z1 = data["z1"]
    z2 = data["z2"]

    data[self.name] = a_int / (z1 + z2)

    ddenominator = - a_int / (z1 + z2)**2
    der[self.name]["vf1"] = der["a_int"]["vf1"] / (z1 + z2) + ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    der[self.name]["arho1"] = ddenominator * der["z1"]["arho1"]
    der[self.name]["arhou1"] = ddenominator * der["z1"]["arhou1"]
    der[self.name]["arhoE1"] = ddenominator * der["z1"]["arhoE1"]
    der[self.name]["arho2"] = ddenominator * der["z2"]["arho2"]
    der[self.name]["arhou2"] = ddenominator * der["z2"]["arhou2"]
    der[self.name]["arhoE2"] = ddenominator * der["z2"]["arhoE2"]
