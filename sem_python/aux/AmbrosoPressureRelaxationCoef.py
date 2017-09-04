from .AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

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

    vf2 = 1 - vf1
    dvf2_dvf1 = 0 * vf2 - 1

    denominator = self.pressure_relaxation_time * (p1 + p2)
    data[self.name] = vf1 * vf2 / denominator
    dp_relax_ddenominator = - vf1 * vf2 / denominator / denominator
    pdp_relax_pdvf1 = (vf2 + vf1 * dvf2_dvf1) / denominator
    dp_relax_dp1 = dp_relax_ddenominator * self.pressure_relaxation_time
    dp_relax_dp2 = dp_relax_ddenominator * self.pressure_relaxation_time

    der[self.name]["vf1"] = pdp_relax_pdvf1 + dp_relax_dp1 * der["p1"]["vf1"] + dp_relax_dp2 * der["p2"]["vf1"]
    der[self.name]["arho1"] = dp_relax_dp1 * der["p1"]["arho1"]
    der[self.name]["arho2"] = dp_relax_dp2 * der["p2"]["arho2"]
    der[self.name]["arhou1"] = dp_relax_dp1 * der["p1"]["arhou1"]
    der[self.name]["arhou2"] = dp_relax_dp2 * der["p2"]["arhou2"]
    der[self.name]["arhoE1"] = dp_relax_dp1 * der["p1"]["arhoE1"]
    der[self.name]["arhoE2"] = dp_relax_dp2 * der["p2"]["arhoE2"]
