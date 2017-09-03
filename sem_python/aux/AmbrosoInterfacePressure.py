from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

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
    der[self.name]["arho1"] = dpI_dp1 * der["p1"]["arho1"] + dpI_du_relax * der["u_relax"]["arho1"]
    der[self.name]["arhou1"] = dpI_dp1 * der["p1"]["arhou1"] + dpI_du_relax * der["u_relax"]["arhou1"]
    der[self.name]["arhoE1"] = dpI_dp1 * der["p1"]["arhoE1"] + dpI_du_relax * der["u_relax"]["arhoE1"]
    der[self.name]["arho2"] = dpI_dp2 * der["p2"]["arho2"] + dpI_du_relax * der["u_relax"]["arho2"]
    der[self.name]["arhou2"] = dpI_dp2 * der["p2"]["arhou2"] + dpI_du_relax * der["u_relax"]["arhou2"]
    der[self.name]["arhoE2"] = dpI_dp2 * der["p2"]["arhoE2"] + dpI_du_relax * der["u_relax"]["arhoE2"]
