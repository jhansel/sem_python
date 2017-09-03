import numpy as np

from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

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
    self.name = self.viscflux_vf

  def compute(self, data, der):
    l = data["visccoef_vf"] * data[self.grad_vf]
    dl_dvisccoef = data[self.grad_vf]
    dl_dgrad_vf = data["visccoef_vf"]

    dl_dvf1 = dl_dvisccoef * der["visccoef_vf"]["vf1"]
    dl_darhoA1 = dl_dvisccoef * der["visccoef_vf"]["arhoA1"]
    dl_darhouA1 = dl_dvisccoef * der["visccoef_vf"]["arhouA1"]
    dl_darhoEA1 = dl_dvisccoef * der["visccoef_vf"]["arhoEA1"]
    dl_darhoA2 = dl_dvisccoef * der["visccoef_vf"]["arhoA2"]
    dl_darhouA2 = dl_dvisccoef * der["visccoef_vf"]["arhouA2"]
    dl_darhoEA2 = dl_dvisccoef * der["visccoef_vf"]["arhoEA2"]
    dl_dgrad_vf1 = dl_dgrad_vf * self.dgrad_vf_dgrad_vf1

    data[self.name] = l
    der[self.name]["vf1"] = dl_dvf1
    der[self.name]["arhoA1"] = dl_darhoA1
    der[self.name]["arhouA1"] = dl_darhouA1
    der[self.name]["arhoEA1"] = dl_darhoEA1
    der[self.name]["arhoA2"] = dl_darhoA2
    der[self.name]["arhouA2"] = dl_darhouA2
    der[self.name]["arhoEA2"] = dl_darhoEA2
    der[self.name]["grad_vf1"] = dl_dgrad_vf1
