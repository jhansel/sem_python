from numpy import sign

from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfacePressureParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryInterfacePressure(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "pI"

  def compute(self, data, der):
    pI_bar = data["pI_bar"]
    grad_vf1 = data["grad_vf1"]
    u1 = data["u1"]
    u2 = data["u2"]
    z1 = data["z1"]
    z2 = data["z2"]

    numerator = sign(grad_vf1) * z1 * z2 * (u2 - u1)
    denominator = z1 + z2
    data[self.name] = pI_bar + numerator / denominator

    dnumerator = 1.0 / denominator
    ddenominator = -numerator / denominator**2
    der[self.name]["vf1"] = der["pI_bar"]["vf1"] + dnumerator * sign(grad_vf1) * (der["z1"]["vf1"] * z2 + z1 * der["z2"]["vf1"]) * (u2 - u1) \
      + ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    der[self.name]["arhoA1"] = der["pI_bar"]["arhoA1"] + dnumerator * sign(grad_vf1) * z2 * (der["z1"]["arhoA1"] * (u2 - u1) - z1 * der["u1"]["arhoA1"]) \
      + ddenominator * der["z1"]["arhoA1"]
    der[self.name]["arhouA1"] = der["pI_bar"]["arhouA1"] + dnumerator * sign(grad_vf1) * z2 * (der["z1"]["arhouA1"] * (u2 - u1) - z1 * der["u1"]["arhouA1"]) \
      + ddenominator * der["z1"]["arhouA1"]
    der[self.name]["arhoEA1"] = der["pI_bar"]["arhoEA1"] + dnumerator * sign(grad_vf1) * der["z1"]["arhoEA1"] * z2 * (u2 - u1) \
      + ddenominator * der["z1"]["arhoEA1"]
    der[self.name]["arhoA2"] = der["pI_bar"]["arhoA2"] + dnumerator * sign(grad_vf1) * z1 * (der["z2"]["arhoA2"] * (u2 - u1) + z2 * der["u2"]["arhoA2"]) \
      + ddenominator * der["z2"]["arhoA2"]
    der[self.name]["arhouA2"] = der["pI_bar"]["arhouA2"] + dnumerator * sign(grad_vf1) * z1 * (der["z2"]["arhouA2"] * (u2 - u1) + z2 * der["u2"]["arhouA2"]) \
      + ddenominator * der["z2"]["arhouA2"]
    der[self.name]["arhoEA2"] = der["pI_bar"]["arhoEA2"] + dnumerator * sign(grad_vf1) * z1 * der["z2"]["arhoEA2"] * (u2 - u1) \
      + ddenominator * der["z2"]["arhoEA2"]
