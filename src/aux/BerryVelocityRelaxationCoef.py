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

    dvf1 = third * (der["p_relax"]["vf1"] * z1 * z2 + p_relax * (der["z1"]["vf1"] * z2 + z1 * der["z2"]["vf1"]))
    darho1 = third * (der["p_relax"]["arho1"] * z1 * z2 + p_relax * der["z1"]["arho1"] * z2)
    darhou1 = third * (der["p_relax"]["arhou1"] * z1 * z2 + p_relax * der["z1"]["arhou1"] * z2)
    darhoE1 = third * (der["p_relax"]["arhoE1"] * z1 * z2 + p_relax * der["z1"]["arhoE1"] * z2)
    darho2 = third * (der["p_relax"]["arho2"] * z1 * z2 + p_relax * z1 * der["z2"]["arho2"])
    darhou2 = third * (der["p_relax"]["arhou2"] * z1 * z2 + p_relax * z1 * der["z2"]["arhou2"])
    darhoE2 = third * (der["p_relax"]["arhoE2"] * z1 * z2 + p_relax * z1 * der["z2"]["arhoE2"])

    der[self.name] = {"vf1": dvf1, "arho1": darho1, "arhou1": darhou1, "arhoE1": darhoE1,
      "arho2": darho2, "arhou2": darhou2, "arhoE2": darhoE2}
