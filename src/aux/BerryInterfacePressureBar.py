from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfacePressureBarParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryInterfacePressureBar(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    p1 = data["p1"]
    p2 = data["p2"]
    z1 = data["z1"]
    z2 = data["z2"]

    numerator = (z1 * p2 + z2 * p1)
    denominator = z1 + z2
    data["pI_bar"] = numerator / denominator

    dnumerator = 1.0 / denominator
    ddenominator = -numerator / denominator**2
    dvf1 = dnumerator * (der["z1"]["vf1"] * p2 + z1 * der["p2"]["vf1"] + der["z2"]["vf1"] * p1 + z2 * der["p1"]["vf1"]) \
      + ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    darho1 = dnumerator * (der["z1"]["arho1"] * p2 + z2 * der["p1"]["arho1"]) \
      + ddenominator * der["z1"]["arho1"]
    darhou1 = dnumerator * (der["z1"]["arhou1"] * p2 + z2 * der["p1"]["arhou1"]) \
      + ddenominator * der["z1"]["arhou1"]
    darhoE1 = dnumerator * (der["z1"]["arhoE1"] * p2 + z2 * der["p1"]["arhoE1"]) \
      + ddenominator * der["z1"]["arhoE1"]
    darho2 = dnumerator * (z1 * der["p2"]["arho2"] + der["z2"]["arho2"] * p1) \
      + ddenominator * der["z2"]["arho2"]
    darhou2 = dnumerator * (z1 * der["p2"]["arhou2"] + der["z2"]["arhou2"] * p1) \
      + ddenominator * der["z2"]["arhou2"]
    darhoE2 = dnumerator * (z1 * der["p2"]["arhoE2"] + der["z2"]["arhoE2"] * p1) \
      + ddenominator * der["z2"]["arhoE2"]

    der["pI_bar"] = {"vf1": dvf1, "arho1": darho1, "arhou1": darhou1, "arhoE1": darhoE1,
      "arho2": darho2, "arhou2": darhou2, "arhoE2": darhoE2}
