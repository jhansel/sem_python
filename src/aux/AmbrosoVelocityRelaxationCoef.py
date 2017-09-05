from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

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

    der[self.name]["aA1"] = du_relax_dT1 * der["T1"]["aA1"] + du_relax_dT2 * der["T2"]["aA1"]
    der[self.name]["arhoA1"] = du_relax_dT1 * der["T1"]["arhoA1"] + du_relax_dbeta * der["beta"]["arhoA1"]
    der[self.name]["arhoA2"] = du_relax_dT2 * der["T2"]["arhoA2"] + du_relax_dbeta * der["beta"]["arhoA2"]
    der[self.name]["arhouA1"] = du_relax_dT1 * der["T1"]["arhouA1"]
    der[self.name]["arhouA2"] = du_relax_dT2 * der["T2"]["arhouA2"]
    der[self.name]["arhoEA1"] = du_relax_dT1 * der["T1"]["arhoEA1"]
    der[self.name]["arhoEA2"] = du_relax_dT2 * der["T2"]["arhoEA2"]
