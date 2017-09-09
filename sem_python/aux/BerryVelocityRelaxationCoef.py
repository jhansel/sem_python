from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryVelocityRelaxationCoefParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class BerryVelocityRelaxationCoef(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "u_relax"

  def compute(self, data, der):
    p_relax = data["p_relax"]
    z1 = data["z1"]
    z2 = data["z2"]

    third = 1.0 / 3.0
    data[self.name] = third * p_relax * z1 * z2

    der[self.name]["aA1"] = third * (der["p_relax"]["aA1"] * z1 * z2 + p_relax * (der["z1"]["aA1"] * z2 + z1 * der["z2"]["aA1"]))
    der[self.name]["arhoA1"] = third * (der["p_relax"]["arhoA1"] * z1 * z2 + p_relax * der["z1"]["arhoA1"] * z2)
    der[self.name]["arhouA1"] = third * (der["p_relax"]["arhouA1"] * z1 * z2 + p_relax * der["z1"]["arhouA1"] * z2)
    der[self.name]["arhoEA1"] = third * (der["p_relax"]["arhoEA1"] * z1 * z2 + p_relax * der["z1"]["arhoEA1"] * z2)
    der[self.name]["arhoA2"] = third * (der["p_relax"]["arhoA2"] * z1 * z2 + p_relax * z1 * der["z2"]["arhoA2"])
    der[self.name]["arhouA2"] = third * (der["p_relax"]["arhouA2"] * z1 * z2 + p_relax * z1 * der["z2"]["arhouA2"])
    der[self.name]["arhoEA2"] = third * (der["p_relax"]["arhoEA2"] * z1 * z2 + p_relax * z1 * der["z2"]["arhoEA2"])
