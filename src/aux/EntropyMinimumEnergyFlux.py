import numpy as np

from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

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
    self.name = self.viscflux_arhoEA

  def compute(self, data, der):
    h = data[self.vf] * data[self.visccoef_arhoEA] * data[self.grad_rhoe] \
      + data[self.rhoe] * data[self.viscflux_vf] \
      - 0.5 * data[self.u] ** 2 * data[self.viscflux_arhoA] \
      + data[self.u] * data[self.viscflux_arhouA]

    dh_dvf1 = (der[self.vf]["vf1"] * data[self.visccoef_arhoEA] \
      + data[self.vf] * der[self.visccoef_arhoEA]["vf1"]) * data[self.grad_rhoe] \
      + der[self.rhoe]["vf1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["vf1"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["vf1"] \
      + data[self.u] * der[self.viscflux_arhouA]["vf1"]
    dh_darhoA1 = data[self.vf] * der[self.visccoef_arhoEA]["arhoA1"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhoA1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhoA1"] \
      - data[self.u] * der[self.u]["arhoA1"] * data[self.viscflux_arhoA] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["arhoA1"] \
      + der[self.u]["arhoA1"] * data[self.viscflux_arhouA] + data[self.u] * der[self.viscflux_arhouA]["arhoA1"]
    dh_darhouA1 = data[self.vf] * der[self.visccoef_arhoEA]["arhouA1"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhouA1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhouA1"] \
      - data[self.u] * der[self.u]["arhouA1"] * data[self.viscflux_arhoA] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["arhouA1"] \
      + der[self.u]["arhouA1"] * data[self.viscflux_arhouA] + data[self.u] * der[self.viscflux_arhouA]["arhouA1"]
    dh_darhoEA1 = data[self.vf] * der[self.visccoef_arhoEA]["arhoEA1"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhoEA1"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhoEA1"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["arhoEA1"] \
      + data[self.u] * der[self.viscflux_arhouA]["arhoEA1"]
    dh_darhoA2 = data[self.vf] * der[self.visccoef_arhoEA]["arhoA2"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhoA2"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhoA2"] \
      - data[self.u] * der[self.u]["arhoA2"] * data[self.viscflux_arhoA] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["arhoA2"] \
      + der[self.u]["arhoA2"] * data[self.viscflux_arhouA] + data[self.u] * der[self.viscflux_arhouA]["arhoA2"]
    dh_darhouA2 = data[self.vf] * der[self.visccoef_arhoEA]["arhouA2"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhouA2"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhouA2"] \
      - data[self.u] * der[self.u]["arhouA2"] * data[self.viscflux_arhoA] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["arhouA2"] \
      + der[self.u]["arhouA2"] * data[self.viscflux_arhouA] + data[self.u] * der[self.viscflux_arhouA]["arhouA2"]
    dh_darhoEA2 = data[self.vf] * der[self.visccoef_arhoEA]["arhoEA2"] * data[self.grad_rhoe] \
      + der[self.rhoe]["arhoEA2"] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf]["arhoEA2"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["arhoEA2"] \
      + data[self.u] * der[self.viscflux_arhouA]["arhoEA2"]
    dh_dgrad_vf1 = data[self.vf] * data[self.visccoef_arhoEA] * der[self.grad_rhoe]["grad_vf1"] \
      + data[self.rhoe] * der[self.viscflux_vf]["grad_vf1"] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA]["grad_vf1"] \
      + data[self.u] * der[self.viscflux_arhouA]["grad_vf1"]
    dh_dgrad_arhoA = data[self.vf] * data[self.visccoef_arhoEA] * der[self.grad_rhoe][self.grad_arhoA] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arhoA][self.grad_arhoA] \
      + data[self.u] * der[self.viscflux_arhouA][self.grad_arhoA]
    dh_dgrad_arhouA = data[self.vf] * data[self.visccoef_arhoEA] * der[self.grad_rhoe][self.grad_arhouA] \
      + data[self.u] * der[self.viscflux_arhouA][self.grad_arhouA]
    dh_dgrad_arhoEA = data[self.vf] * data[self.visccoef_arhoEA] * der[self.grad_rhoe][self.grad_arhoEA]

    data[self.name] = h
    der[self.name]["vf1"] = dh_dvf1
    der[self.name]["arhoA1"] = dh_darhoA1
    der[self.name]["arhouA1"] = dh_darhouA1
    der[self.name]["arhoEA1"] = dh_darhoEA1
    der[self.name]["arhoA2"] = dh_darhoA2
    der[self.name]["arhouA2"] = dh_darhouA2
    der[self.name]["arhoEA2"] = dh_darhoEA2
    der[self.name]["grad_vf1"] = dh_dgrad_vf1
    der[self.name][self.grad_arhoA] = dh_dgrad_arhoA
    der[self.name][self.grad_arhouA] = dh_dgrad_arhouA
    der[self.name][self.grad_arhoEA] = dh_dgrad_arhoEA
