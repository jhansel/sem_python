from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoPressureRelaxationCoefParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerFloatParameter("pressure_relaxation_time", "Relaxation time for pressures")

class AmbrosoPressureRelaxationCoef(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "p_relax"
    self.pressure_relaxation_time = params.get("pressure_relaxation_time")

  def compute(self, data, der):
    vf1 = data["vf1"]
    p1 = data["p1"]
    p2 = data["p2"]

    denominator = self.pressure_relaxation_time * (p1 + p2)
    data[self.name] = vf1 * (1 - vf1) / denominator
    dp_relax_ddenominator = - vf1 * (1 - vf1) / denominator / denominator
    pdp_relax_pdaA1 = (1 - 2 * vf1) * der["vf1"]["aA1"] / denominator
    dp_relax_dp1 = dp_relax_ddenominator * self.pressure_relaxation_time
    dp_relax_dp2 = dp_relax_ddenominator * self.pressure_relaxation_time

    der[self.name]["aA1"] = pdp_relax_pdaA1 + dp_relax_dp1 * der["p1"]["aA1"] + dp_relax_dp2 * der["p2"]["aA1"]
    der[self.name]["arhoA1"] = dp_relax_dp1 * der["p1"]["arhoA1"]
    der[self.name]["arhoA2"] = dp_relax_dp2 * der["p2"]["arhoA2"]
    der[self.name]["arhouA1"] = dp_relax_dp1 * der["p1"]["arhouA1"]
    der[self.name]["arhouA2"] = dp_relax_dp2 * der["p2"]["arhouA2"]
    der[self.name]["arhoEA1"] = dp_relax_dp1 * der["p1"]["arhoEA1"]
    der[self.name]["arhoEA2"] = dp_relax_dp2 * der["p2"]["arhoEA2"]
