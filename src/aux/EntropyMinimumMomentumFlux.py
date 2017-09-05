import numpy as np

from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class EntropyMinimumMomentumFluxParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

## Entropy minimum viscous flux for the momentum equation
#
# For phase \f$k\f$, this flux is computed as
# \f[
#   g_k \equiv \alpha_k\nu_k\rho_k\nabla u_k + u_k f_k ,
# \f]
# where \f$f_k\f$ is the viscous flux for the mass equation.
#
class EntropyMinimumMomentumFlux(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.viscflux_arhouA

  def compute(self, data, der):
    g = data[self.vf] * data[self.visccoef_arhouA] * data[self.rho] * data["A"] * data[self.grad_u] \
      + data[self.u] * data[self.viscflux_arhoA]

    dg_daA1 = (der[self.vf]["aA1"] * data[self.visccoef_arhouA] * data[self.rho] \
      + data[self.vf] * der[self.visccoef_arhouA]["aA1"] * data[self.rho] \
      + data[self.vf] * data[self.visccoef_arhouA] * der[self.rho]["aA1"]) * data["A"] * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arhoA]["aA1"]
    dg_darhoA1 = data[self.vf] * (der[self.visccoef_arhouA]["arhoA1"] * data[self.rho] \
      + data[self.visccoef_arhouA] * der[self.rho]["arhoA1"]) * data["A"] * data[self.grad_u] \
      + der[self.u]["arhoA1"] * data[self.viscflux_arhoA] + data[self.u] * der[self.viscflux_arhoA]["arhoA1"]
    dg_darhouA1 = data[self.vf] * der[self.visccoef_arhouA]["arhouA1"] * data[self.rho] * data["A"] * data[self.grad_u] \
      + der[self.u]["arhouA1"] * data[self.viscflux_arhoA] + data[self.u] * der[self.viscflux_arhoA]["arhouA1"]
    dg_darhoEA1 = data[self.vf] * der[self.visccoef_arhouA]["arhoEA1"] * data[self.rho] * data["A"] * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arhoA]["arhoEA1"]
    dg_darhoA2 = data[self.vf] * (der[self.visccoef_arhouA]["arhoA2"] * data[self.rho] \
      + data[self.visccoef_arhouA] * der[self.rho]["arhoA2"]) * data["A"] * data[self.grad_u] \
      + der[self.u]["arhoA2"] * data[self.viscflux_arhoA] + data[self.u] * der[self.viscflux_arhoA]["arhoA2"]
    dg_darhouA2 = data[self.vf] * der[self.visccoef_arhouA]["arhouA2"] * data[self.rho] * data["A"] * data[self.grad_u] \
      + der[self.u]["arhouA2"] * data[self.viscflux_arhoA] + data[self.u] * der[self.viscflux_arhoA]["arhouA2"]
    dg_darhoEA2 = data[self.vf] * der[self.visccoef_arhouA]["arhoEA2"] * data[self.rho] * data["A"] * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arhoA]["arhoEA2"]
    dg_dgrad_aA1 = data[self.u] * der[self.viscflux_arhoA]["grad_aA1"]
    dg_dgrad_arhoA = data[self.vf] * data[self.visccoef_arhouA] * data[self.rho] * data["A"] * der[self.grad_u][self.grad_arhoA] \
      + data[self.u] * der[self.viscflux_arhoA][self.grad_arhoA]
    dg_dgrad_arhouA = data[self.vf] * data[self.visccoef_arhouA] * data[self.rho] * data["A"] * der[self.grad_u][self.grad_arhouA]

    data[self.name] = g
    der[self.name]["aA1"] = dg_daA1
    der[self.name]["arhoA1"] = dg_darhoA1
    der[self.name]["arhouA1"] = dg_darhouA1
    der[self.name]["arhoEA1"] = dg_darhoEA1
    der[self.name]["arhoA2"] = dg_darhoA2
    der[self.name]["arhouA2"] = dg_darhouA2
    der[self.name]["arhoEA2"] = dg_darhoEA2
    der[self.name]["grad_aA1"] = dg_dgrad_aA1
    der[self.name][self.grad_arhoA] = dg_dgrad_arhoA
    der[self.name][self.grad_arhouA] = dg_dgrad_arhouA
