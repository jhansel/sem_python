from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoVelocityRelaxationCoefParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)

class AmbrosoVelocityRelaxationCoef(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)

  def compute(self, data, der):
    beta = data["beta"]
    T1 = data["T1"]
    T2 = data["T2"]

    denominator = beta * T1 + (1 - beta) * T2
    data["u_relax"] = (1 - beta) * T2 / denominator
    du_relax_dT1 = - (1 - beta) * T2 / denominator / denominator * beta
    du_relax_dT2 = (1 - beta) / denominator - (1 - beta) * T2 / denominator / denominator * (1 - beta)
    du_relax_dbeta = - T2 / denominator - (1 - beta) * T2 / denominator / denominator * (T1 - T2)

    du_relax_dvf1 = du_relax_dT1 * der["T1"]["vf1"] + du_relax_dT2 * der["T2"]["vf1"]
    du_relax_darho1 = du_relax_dT1 * der["T1"]["arho1"] + du_relax_dbeta * der["beta"]["arho1"]
    du_relax_darho2 = du_relax_dT2 * der["T2"]["arho2"] + du_relax_dbeta * der["beta"]["arho2"]
    du_relax_darhou1 = du_relax_dT1 * der["T1"]["arhou1"]
    du_relax_darhou2 = du_relax_dT2 * der["T2"]["arhou2"]
    du_relax_darhoE1 = du_relax_dT1 * der["T1"]["arhoE1"]
    du_relax_darhoE2 = du_relax_dT2 * der["T2"]["arhoE2"]

    der["u_relax"] = {"vf1": du_relax_dvf1,
                 "arho1": du_relax_darho1, "arhou1": du_relax_darhou1, "arhoE1": du_relax_darhoE1,
                 "arho2": du_relax_darho2, "arhou2": du_relax_darhou2, "arhoE2": du_relax_darhoE2}
