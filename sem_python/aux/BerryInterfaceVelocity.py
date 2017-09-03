from numpy import sign

from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfaceVelocityParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryInterfaceVelocity(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "uI"

  def compute(self, data, der):
    uI_bar = data["uI_bar"]
    grad_vf1 = data["grad_vf1"]
    p1 = data["p1"]
    p2 = data["p2"]
    z1 = data["z1"]
    z2 = data["z2"]

    numerator = sign(grad_vf1) * (p2 - p1)
    denominator = z1 + z2
    data[self.name] = uI_bar + numerator / denominator

    dnumerator = 1.0 / denominator
    ddenominator = -numerator / denominator**2
    der[self.name]["vf1"] = der["uI_bar"]["vf1"] + dnumerator * sign(grad_vf1) * (der["p2"]["vf1"] - der["p1"]["vf1"]) \
      + ddenominator * (der["z1"]["vf1"] + der["z2"]["vf1"])
    der[self.name]["arho1"] = der["uI_bar"]["arho1"] - dnumerator * sign(grad_vf1) * der["p1"]["arho1"] + ddenominator * der["z1"]["arho1"]
    der[self.name]["arhou1"] = der["uI_bar"]["arhou1"] - dnumerator * sign(grad_vf1) * der["p1"]["arhou1"] + ddenominator * der["z1"]["arhou1"]
    der[self.name]["arhoE1"] = der["uI_bar"]["arhoE1"] - dnumerator * sign(grad_vf1) * der["p1"]["arhoE1"] + ddenominator * der["z1"]["arhoE1"]
    der[self.name]["arho2"] = der["uI_bar"]["arho2"] + dnumerator * sign(grad_vf1) * der["p2"]["arho2"] + ddenominator * der["z2"]["arho2"]
    der[self.name]["arhou2"] = der["uI_bar"]["arhou2"] + dnumerator * sign(grad_vf1) * der["p2"]["arhou2"] + ddenominator * der["z2"]["arhou2"]
    der[self.name]["arhoE2"] = der["uI_bar"]["arhoE2"] + dnumerator * sign(grad_vf1) * der["p2"]["arhoE2"] + ddenominator * der["z2"]["arhoE2"]
