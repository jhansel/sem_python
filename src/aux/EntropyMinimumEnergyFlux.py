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
    dh_darho = data[self.vf] * der[self.visccoef_arhoE][self.arho] * data[self.grad_rhoe] \
      + der[self.rhoe][self.arho] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf][self.arho] \
      - data[self.u] * der[self.u][self.arho] * data[self.viscflux_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho][self.arho] \
      + der[self.u][self.arho] * data[self.viscflux_arhou] + data[self.u] * der[self.viscflux_arhou][self.arho]
    dh_darhou = data[self.vf] * der[self.visccoef_arhoE][self.arhou] * data[self.grad_rhoe] \
      + der[self.rhoe][self.arhou] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf][self.arhou] \
      - data[self.u] * der[self.u][self.arhou] * data[self.viscflux_arho] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho][self.arhou] \
      + der[self.u][self.arhou] * data[self.viscflux_arhou] + data[self.u] * der[self.viscflux_arhou][self.arhou]
    dh_darhoE = data[self.vf] * der[self.visccoef_arhoE][self.arhoE] * data[self.grad_rhoe] \
      + der[self.rhoe][self.arhoE] * data[self.viscflux_vf] + data[self.rhoe] * der[self.viscflux_vf][self.arhoE] \
      - 0.5 * data[self.u] ** 2 * der[self.viscflux_arho][self.arhoE] \
      + data[self.u] * der[self.viscflux_arhou][self.arhoE]
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
    der[self.name] = {"vf1" : dh_dvf1, self.arho : dh_darho, self.arhou : dh_darhou,
      self.arhoE : dh_darhoE, "grad_vf1" : dh_dgrad_vf1, self.grad_arho : dh_dgrad_arho,
      self.grad_arhou : dh_dgrad_arhou, self.grad_arhoE : dh_dgrad_arhoE}
