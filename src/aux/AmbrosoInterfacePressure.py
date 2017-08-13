from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoInterfacePressureParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoInterfacePressure(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    u_relax = data["u_relax"]
    p1 = data["p1"]
    p2 = data["p2"]

    data["pI"] = u_relax * p1 + (1 - u_relax) * p2
    dpI_dp1 = u_relax
    dpI_dp2 = (1 - u_relax)
    dpI_du_relax = p1 - p2

    dpI_dvf1 = dpI_dp1 * der["p1"]["vf1"] + dpI_dp2 * der["p2"]["vf1"] + dpI_du_relax * der["u_relax"]["vf1"]
    dpI_darho1 = dpI_dp1 * der["p1"]["arho1"] + dpI_du_relax * der["u_relax"]["arho1"]
    dpI_darho2 = dpI_dp2 * der["p2"]["arho2"] + dpI_du_relax * der["u_relax"]["arho2"]
    dpI_darhou1 = dpI_dp1 * der["p1"]["arhou1"] + dpI_du_relax * der["u_relax"]["arhou1"]
    dpI_darhoE1 = dpI_dp1 * der["p1"]["arhoE1"] + dpI_du_relax * der["u_relax"]["arhoE1"]
    dpI_darhou2 = dpI_dp2 * der["p2"]["arhou2"] + dpI_du_relax * der["u_relax"]["arhou2"]
    dpI_darhoE2 = dpI_dp2 * der["p2"]["arhoE2"] + dpI_du_relax * der["u_relax"]["arhoE2"]
    der["pI"] = {"vf1": dpI_dvf1,
                 "arho1": dpI_darho1, "arhou1": dpI_darhou1, "arhoE1": dpI_darhoE1,
                 "arho2": dpI_darho2, "arhou2": dpI_darhou2, "arhoE2": dpI_darhoE2}
