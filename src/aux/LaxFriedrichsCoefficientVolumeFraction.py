import numpy as np

from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class LaxFriedrichsCoefficientVolumeFractionParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerFloatParameter("mult", "Multiplier factor for the viscous coefficient", 1.0)

## Lax-Friedrichs artificial viscosity coefficient for the volume fraction equation
#
# This coefficient is computed as
# \f[
#   b_k \equiv \frac{1}{2}\Delta x \left|u_I\right| ,
# \f]
# where \f$\Delta x\f$ is the element width.
#
class LaxFriedrichsCoefficientVolumeFraction(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.name = "visccoef_vf"
    self.mult = params.get("mult")

  def compute(self, data, der):
    b = self.mult * 0.5 * data["dx"] * np.abs(data["uI"])
    db_duI = self.mult * 0.5 * data["dx"] * np.sign(data["uI"])

    db_dvf1 = db_duI * der["uI"]["vf1"]
    db_darhoA1 = db_duI * der["uI"]["arhoA1"]
    db_darhouA1 = db_duI * der["uI"]["arhouA1"]
    db_darhoEA1 = db_duI * der["uI"]["arhoEA1"]
    db_darhoA2 = db_duI * der["uI"]["arhoA2"]
    db_darhouA2 = db_duI * der["uI"]["arhouA2"]
    db_darhoEA2 = db_duI * der["uI"]["arhoEA2"]

    data[self.name] = b
    der[self.name]["vf1"] = db_dvf1
    der[self.name]["arhoA1"] = db_darhoA1
    der[self.name]["arhouA1"] = db_darhouA1
    der[self.name]["arhoEA1"] = db_darhoEA1
    der[self.name]["arhoA2"] = db_darhoA2
    der[self.name]["arhouA2"] = db_darhouA2
    der[self.name]["arhoEA2"] = db_darhoEA2
