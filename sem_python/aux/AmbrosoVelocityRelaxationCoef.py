from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoVelocityRelaxationCoefParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoVelocityRelaxationCoef(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "u_relax"

  def compute(self, data, der):
    beta = data["beta"]
    T1 = data["T1"]
    T2 = data["T2"]

    denominator = beta * T1 + (1 - beta) * T2
    data[self.name] = (1 - beta) * T2 / denominator
    du_relax_dT1 = - (1 - beta) * T2 / denominator / denominator * beta
    du_relax_dT2 = (1 - beta) / denominator - (1 - beta) * T2 / denominator / denominator * (1 - beta)
    du_relax_dbeta = - T2 / denominator - (1 - beta) * T2 / denominator / denominator * (T1 - T2)

    der[self.name]["vf1"] = du_relax_dT1 * der["T1"]["vf1"] + du_relax_dT2 * der["T2"]["vf1"]
    der[self.name]["arho1"] = du_relax_dT1 * der["T1"]["arho1"] + du_relax_dbeta * der["beta"]["arho1"]
    der[self.name]["arho2"] = du_relax_dT2 * der["T2"]["arho2"] + du_relax_dbeta * der["beta"]["arho2"]
    der[self.name]["arhou1"] = du_relax_dT1 * der["T1"]["arhou1"]
    der[self.name]["arhou2"] = du_relax_dT2 * der["T2"]["arhou2"]
    der[self.name]["arhoE1"] = du_relax_dT1 * der["T1"]["arhoE1"]
    der[self.name]["arhoE2"] = du_relax_dT2 * der["T2"]["arhoE2"]
