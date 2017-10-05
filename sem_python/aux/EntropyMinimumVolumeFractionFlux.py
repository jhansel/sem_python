import numpy as np

from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class EntropyMinimumVolumeFractionFluxParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

## Entropy minimum viscous flux for the volume fraction equation
#
# For phase \f$k\f$, this flux is computed as
# \f[
#   l_k \equiv \nu_k \nabla \alpha_k ,
# \f]
# where \f$\nu_k\f$ is the viscous coefficient.
#
class EntropyMinimumVolumeFractionFlux(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.viscflux_aA

  def compute(self, data, der):
    l = data["visccoef_aA1"] * data["A"] * data[self.grad_vf]
    dl_dvisccoef = data["A"] * data[self.grad_vf]
    dl_dgrad_vf = data["visccoef_aA1"] * data["A"]

    dl_daA1 = dl_dvisccoef * der["visccoef_aA1"]["aA1"]
    dl_darhoA1 = dl_dvisccoef * der["visccoef_aA1"]["arhoA1"]
    dl_darhouA1 = dl_dvisccoef * der["visccoef_aA1"]["arhouA1"]
    dl_darhoEA1 = dl_dvisccoef * der["visccoef_aA1"]["arhoEA1"]
    dl_darhoA2 = dl_dvisccoef * der["visccoef_aA1"]["arhoA2"]
    dl_darhouA2 = dl_dvisccoef * der["visccoef_aA1"]["arhouA2"]
    dl_darhoEA2 = dl_dvisccoef * der["visccoef_aA1"]["arhoEA2"]
    dl_dgrad_aA1 = dl_dgrad_vf * der[self.grad_vf]["grad_aA1"]

    data[self.name] = l
    der[self.name]["aA1"] = dl_daA1
    der[self.name]["arhoA1"] = dl_darhoA1
    der[self.name]["arhouA1"] = dl_darhouA1
    der[self.name]["arhoEA1"] = dl_darhoEA1
    der[self.name]["arhoA2"] = dl_darhoA2
    der[self.name]["arhouA2"] = dl_darhouA2
    der[self.name]["arhoEA2"] = dl_darhoEA2
    der[self.name]["grad_aA1"] = dl_dgrad_aA1
