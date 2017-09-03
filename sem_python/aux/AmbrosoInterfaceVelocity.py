from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

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

    der[self.name]["arho1"] = duI_du1 * der["u1"]["arho1"] + duI_dbeta * der["beta"]["arho1"]
    der[self.name]["arho2"] = duI_du2 * der["u2"]["arho2"] + duI_dbeta * der["beta"]["arho2"]
    der[self.name]["arhou1"] = duI_du1 * der["u1"]["arhou1"]
    der[self.name]["arhou2"] = duI_du2 * der["u2"]["arhou2"]
