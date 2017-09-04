import numpy as np

from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class EntropyMinimumEnergyFluxParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

## Entropy minimum viscous flux for the energy equation
#
# For phase \f$k\f$, this flux is computed as
# \f[
#   h_k \equiv \alpha_k\nu_k\nabla(\rho e)_k + (\rho e)_k l_k - \frac{1}{2}u_k^2 f_k + u_k g_k ,
# \f]
# where \f$l_k\f$ is the viscous flux for the volume fraction equation,
# \f$f_k\f$ is the viscous flux for the mass equation,
# \f$g_k\f$ is the viscous flux for the momentum equation,
# and \f$\nu_k\f$ is the viscous coefficient.
#
class EntropyMinimumEnergyFlux(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.viscflux_arhoE

  def compute(self, data, der):
    h = data[self.vf] * data[self.visccoef_arhoE] * data[self.grad_rhoe] \
      + data[self.rhoe] * data[self.viscflux_vf] \
      - 0.5 * data[self.u] ** 2 * data[self.viscflux_arho] \
      + data[self.u] * data[self.viscflux_arhou]

    dh_dvf1 = (der[self.vf]["vf1"] * data[self.visccoef_arhoE] \
      + data[self.vf] * der[self.visccoef_arhoE]["vf1"]) * data[self.grad_rhoe] \
      + der[self.rhoe]["vf1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["vf1"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["vf1"] \
      + data[self.u] * der[self.viscflux_arhou]["vf1"]
    dh_darho1 = data[self.vf] * der[self.visccoef_arhoE]["arho1"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arho1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arho1"] \
      - data[self.u] * der[self.u]["arho1"] * data[self.viscflux_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["arho1"] \
      + der[self.u]["arho1"] * data[self.viscflux_arhou] + data[self.u] * der[self.viscflux_arhou]["arho1"]
    dh_darhou1 = data[self.vf] * der[self.visccoef_arhoE]["arhou1"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhou1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhou1"] \
      - data[self.u] * der[self.u]["arhou1"] * data[self.viscflux_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["arhou1"] \
      + der[self.u]["arhou1"] * data[self.viscflux_arhou] + data[self.u] * der[self.viscflux_arhou]["arhou1"]
    dh_darhoE1 = data[self.vf] * der[self.visccoef_arhoE]["arhoE1"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhoE1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhoE1"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["arhoE1"] \
      + data[self.u] * der[self.viscflux_arhou]["arhoE1"]
    dh_darho2 = data[self.vf] * der[self.visccoef_arhoE]["arho2"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arho2"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arho2"] \
      - data[self.u] * der[self.u]["arho2"] * data[self.viscflux_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["arho2"] \
      + der[self.u]["arho2"] * data[self.viscflux_arhou] + data[self.u] * der[self.viscflux_arhou]["arho2"]
    dh_darhou2 = data[self.vf] * der[self.visccoef_arhoE]["arhou2"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhou2"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhou2"] \
      - data[self.u] * der[self.u]["arhou2"] * data[self.viscflux_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["arhou2"] \
      + der[self.u]["arhou2"] * data[self.viscflux_arhou] + data[self.u] * der[self.viscflux_arhou]["arhou2"]
    dh_darhoE2 = data[self.vf] * der[self.visccoef_arhoE]["arhoE2"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhoE2"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhoE2"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["arhoE2"] \
      + data[self.u] * der[self.viscflux_arhou]["arhoE2"]
    dh_dgrad_vf1 = data[self.vf] * data[self.visccoef_arhoE] * der[self.grad_rhoe]["grad_vf1"] \
      + data[self.rhoe] * der[self.viscflux_vf]["grad_vf1"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho]["grad_vf1"] \
      + data[self.u] * der[self.viscflux_arhou]["grad_vf1"]
    dh_dgrad_arho = data[self.vf] * data[self.visccoef_arhoE] * der[self.grad_rhoe][self.grad_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho][self.grad_arho] \
      + data[self.u] * der[self.viscflux_arhou][self.grad_arho]
    dh_dgrad_arhou = data[self.vf] * data[self.visccoef_arhoE] * der[self.grad_rhoe][self.grad_arhou] \
      + data[self.u] * der[self.viscflux_arhou][self.grad_arhou]
    dh_dgrad_arhoE = data[self.vf] * data[self.visccoef_arhoE] * der[self.grad_rhoe][self.grad_arhoE]

    data[self.name] = h
    der[self.name]["vf1"] = dh_dvf1
    der[self.name]["arho1"] = dh_darho1
    der[self.name]["arhou1"] = dh_darhou1
    der[self.name]["arhoE1"] = dh_darhoE1
    der[self.name]["arho2"] = dh_darho2
    der[self.name]["arhou2"] = dh_darhou2
    der[self.name]["arhoE2"] = dh_darhoE2
    der[self.name]["grad_vf1"] = dh_dgrad_vf1
    der[self.name][self.grad_arho] = dh_dgrad_arho
    der[self.name][self.grad_arhou] = dh_dgrad_arhou
    der[self.name][self.grad_arhoE] = dh_dgrad_arhoE
