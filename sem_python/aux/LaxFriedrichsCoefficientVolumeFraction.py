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
    db_darho1 = db_duI * der["uI"]["arho1"]
    db_darhou1 = db_duI * der["uI"]["arhou1"]
    db_darhoE1 = db_duI * der["uI"]["arhoE1"]
    db_darho2 = db_duI * der["uI"]["arho2"]
    db_darhou2 = db_duI * der["uI"]["arhou2"]
    db_darhoE2 = db_duI * der["uI"]["arhoE2"]

    data[self.name] = b
    der[self.name]["vf1"] = db_dvf1
    der[self.name]["arho1"] = db_darho1
    der[self.name]["arhou1"] = db_darhou1
    der[self.name]["arhoE1"] = db_darhoE1
    der[self.name]["arho2"] = db_darho2
    der[self.name]["arhou2"] = db_darhou2
    der[self.name]["arhoE2"] = db_darhoE2
