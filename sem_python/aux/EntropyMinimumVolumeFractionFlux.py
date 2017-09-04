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
    self.name = self.viscflux_vf

  def compute(self, data, der):
    l = data["visccoef_vf"] * data[self.grad_vf]
    dl_dvisccoef = data[self.grad_vf]
    dl_dgrad_vf = data["visccoef_vf"]

    dl_dvf1 = dl_dvisccoef * der["visccoef_vf"]["vf1"]
    dl_darho1 = dl_dvisccoef * der["visccoef_vf"]["arho1"]
    dl_darhou1 = dl_dvisccoef * der["visccoef_vf"]["arhou1"]
    dl_darhoE1 = dl_dvisccoef * der["visccoef_vf"]["arhoE1"]
    dl_darho2 = dl_dvisccoef * der["visccoef_vf"]["arho2"]
    dl_darhou2 = dl_dvisccoef * der["visccoef_vf"]["arhou2"]
    dl_darhoE2 = dl_dvisccoef * der["visccoef_vf"]["arhoE2"]
    dl_dgrad_vf1 = dl_dgrad_vf * self.dgrad_vf_dgrad_vf1

    data[self.name] = l
    der[self.name]["vf1"] = dl_dvf1
    der[self.name]["arho1"] = dl_darho1
    der[self.name]["arhou1"] = dl_darhou1
    der[self.name]["arhoE1"] = dl_darhoE1
    der[self.name]["arho2"] = dl_darho2
    der[self.name]["arhou2"] = dl_darhou2
    der[self.name]["arhoE2"] = dl_darhoE2
    der[self.name]["grad_vf1"] = dl_dgrad_vf1
