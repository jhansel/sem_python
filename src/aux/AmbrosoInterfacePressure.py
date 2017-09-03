from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoInterfacePressureParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoInterfacePressure(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "pI"

  def compute(self, data, der):
    u_relax = data["u_relax"]
    p1 = data["p1"]
    p2 = data["p2"]

    data[self.name] = u_relax * p1 + (1 - u_relax) * p2
    dpI_dp1 = u_relax
    dpI_dp2 = (1 - u_relax)
    dpI_du_relax = p1 - p2

    der[self.name]["vf1"] = dpI_dp1 * der["p1"]["vf1"] + dpI_dp2 * der["p2"]["vf1"] + dpI_du_relax * der["u_relax"]["vf1"]
    der[self.name]["arhoA1"] = dpI_dp1 * der["p1"]["arhoA1"] + dpI_du_relax * der["u_relax"]["arhoA1"]
    der[self.name]["arhouA1"] = dpI_dp1 * der["p1"]["arhouA1"] + dpI_du_relax * der["u_relax"]["arhouA1"]
    der[self.name]["arhoEA1"] = dpI_dp1 * der["p1"]["arhoEA1"] + dpI_du_relax * der["u_relax"]["arhoEA1"]
    der[self.name]["arhoA2"] = dpI_dp2 * der["p2"]["arhoA2"] + dpI_du_relax * der["u_relax"]["arhoA2"]
    der[self.name]["arhouA2"] = dpI_dp2 * der["p2"]["arhouA2"] + dpI_du_relax * der["u_relax"]["arhouA2"]
    der[self.name]["arhoEA2"] = dpI_dp2 * der["p2"]["arhoEA2"] + dpI_du_relax * der["u_relax"]["arhoEA2"]
