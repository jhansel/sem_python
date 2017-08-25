from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfaceVelocityBarParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryInterfaceVelocityBar(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "uI_bar"

  def compute(self, data, der):
    u1 = data["u1"]
    u2 = data["u2"]
    z1 = data["z1"]
    z2 = data["z2"]

    numerator = (z1 * u1 + z2 * u2)
    data[self.name] = numerator / (z1 + z2)

    ddenominator = -1.0 / (z1 + z2)**2
    der[self.name]["vf1"] = (der["z1"]["vf1"] * u1 + der["z2"]["vf1"] * u2) / (z1 + z2) \
      + numerator * ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    der[self.name]["arho1"] = (der["z1"]["arho1"] * u1 + z1 * der["u1"]["arho1"]) / (z1 + z2) \
      + numerator * ddenominator * der["z1"]["arho1"]
    der[self.name]["arhou1"] = (der["z1"]["arhou1"] * u1 + z1 * der["u1"]["arhou1"]) / (z1 + z2) \
      + numerator * ddenominator * der["z1"]["arhou1"]
    der[self.name]["arhoE1"] = der["z1"]["arhoE1"] * u1 / (z1 + z2) + numerator * ddenominator * der["z1"]["arhoE1"]
    der[self.name]["arho2"] = (der["z2"]["arho2"] * u2 + z2 * der["u2"]["arho2"]) / (z1 + z2) \
      + numerator * ddenominator * der["z2"]["arho2"]
    der[self.name]["arhou2"] = (der["z2"]["arhou2"] * u2 + z2 * der["u2"]["arhou2"]) / (z1 + z2) \
      + numerator * ddenominator * der["z2"]["arhou2"]
    der[self.name]["arhoE2"] = der["z2"]["arhoE2"] * u2 / (z1 + z2) + numerator * ddenominator * der["z2"]["arhoE2"]
