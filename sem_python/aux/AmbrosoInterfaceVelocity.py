from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoInterfaceVelocityParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoInterfaceVelocity(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "uI"

  def compute(self, data, der):
    beta = data["beta"]
    u1 = data["u1"]
    u2 = data["u2"]

    data[self.name] = beta * u1 + (1 - beta) * u2
    duI_du1 = beta
    duI_du2 = (1 - beta)
    duI_dbeta = u1 - u2

    der[self.name]["arhoA1"] = duI_du1 * der["u1"]["arhoA1"] + duI_dbeta * der["beta"]["arhoA1"]
    der[self.name]["arhoA2"] = duI_du2 * der["u2"]["arhoA2"] + duI_dbeta * der["beta"]["arhoA2"]
    der[self.name]["arhouA1"] = duI_du1 * der["u1"]["arhouA1"]
    der[self.name]["arhouA2"] = duI_du2 * der["u2"]["arhouA2"]
