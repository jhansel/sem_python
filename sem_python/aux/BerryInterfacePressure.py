from numpy import sign

from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

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
    der[self.name]["arho1"] = der["pI_bar"]["arho1"] + dnumerator * sign(grad_vf1) * z2 * (der["z1"]["arho1"] * (u2 - u1) - z1 * der["u1"]["arho1"]) \
      + ddenominator * der["z1"]["arho1"]
    der[self.name]["arhou1"] = der["pI_bar"]["arhou1"] + dnumerator * sign(grad_vf1) * z2 * (der["z1"]["arhou1"] * (u2 - u1) - z1 * der["u1"]["arhou1"]) \
      + ddenominator * der["z1"]["arhou1"]
    der[self.name]["arhoE1"] = der["pI_bar"]["arhoE1"] + dnumerator * sign(grad_vf1) * der["z1"]["arhoE1"] * z2 * (u2 - u1) \
      + ddenominator * der["z1"]["arhoE1"]
    der[self.name]["arho2"] = der["pI_bar"]["arho2"] + dnumerator * sign(grad_vf1) * z1 * (der["z2"]["arho2"] * (u2 - u1) + z2 * der["u2"]["arho2"]) \
      + ddenominator * der["z2"]["arho2"]
    der[self.name]["arhou2"] = der["pI_bar"]["arhou2"] + dnumerator * sign(grad_vf1) * z1 * (der["z2"]["arhou2"] * (u2 - u1) + z2 * der["u2"]["arhou2"]) \
      + ddenominator * der["z2"]["arhou2"]
    der[self.name]["arhoE2"] = der["pI_bar"]["arhoE2"] + dnumerator * sign(grad_vf1) * z1 * der["z2"]["arhoE2"] * (u2 - u1) \
      + ddenominator * der["z2"]["arhoE2"]
