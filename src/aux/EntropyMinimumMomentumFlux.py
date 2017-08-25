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
    self.name = self.viscflux_arhou

  def compute(self, data, der):
    g = data[self.vf] * data[self.visccoef_arhou] * data[self.rho] * data[self.grad_u] \
      + data[self.u] * data[self.viscflux_arho]

    dg_dvf1 = (der[self.vf]["vf1"] * data[self.visccoef_arhou] * data[self.rho] \
      + data[self.vf] * der[self.visccoef_arhou]["vf1"] * data[self.rho] \
      + data[self.vf] * data[self.visccoef_arhou] * der[self.rho]["vf1"]) * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arho]["vf1"]
    dg_darho1 = data[self.vf] * (der[self.visccoef_arhou]["arho1"] * data[self.rho] \
      + data[self.visccoef_arhou] * der[self.rho]["arho1"]) * data[self.grad_u] \
      + der[self.u]["arho1"] * data[self.viscflux_arho] + data[self.u] * der[self.viscflux_arho]["arho1"]
    dg_darhou1 = data[self.vf] * der[self.visccoef_arhou]["arhou1"] * data[self.rho] * data[self.grad_u] \
      + der[self.u]["arhou1"] * data[self.viscflux_arho] + data[self.u] * der[self.viscflux_arho]["arhou1"]
    dg_darhoE1 = data[self.vf] * der[self.visccoef_arhou]["arhoE1"] * data[self.rho] * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arho]["arhoE1"]
    dg_darho2 = data[self.vf] * (der[self.visccoef_arhou]["arho2"] * data[self.rho] \
      + data[self.visccoef_arhou] * der[self.rho]["arho2"]) * data[self.grad_u] \
      + der[self.u]["arho2"] * data[self.viscflux_arho] + data[self.u] * der[self.viscflux_arho]["arho2"]
    dg_darhou2 = data[self.vf] * der[self.visccoef_arhou]["arhou2"] * data[self.rho] * data[self.grad_u] \
      + der[self.u]["arhou2"] * data[self.viscflux_arho] + data[self.u] * der[self.viscflux_arho]["arhou2"]
    dg_darhoE2 = data[self.vf] * der[self.visccoef_arhou]["arhoE2"] * data[self.rho] * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arho]["arhoE2"]
    dg_dgrad_vf1 = data[self.u] * der[self.viscflux_arho]["grad_vf1"]
    dg_dgrad_arho = data[self.vf] * data[self.visccoef_arhou] * data[self.rho] * der[self.grad_u][self.grad_arho] \
      + data[self.u] * der[self.viscflux_arho][self.grad_arho]
    dg_dgrad_arhou = data[self.vf] * data[self.visccoef_arhou] * data[self.rho] * der[self.grad_u][self.grad_arhou]

    data[self.name] = g
    der[self.name]["vf1"] = dg_dvf1
    der[self.name]["arho1"] = dg_darho1
    der[self.name]["arhou1"] = dg_darhou1
    der[self.name]["arhoE1"] = dg_darhoE1
    der[self.name]["arho2"] = dg_darho2
    der[self.name]["arhou2"] = dg_darhou2
    der[self.name]["arhoE2"] = dg_darhoE2
    der[self.name]["grad_vf1"] = dg_dgrad_vf1
    der[self.name][self.grad_arho] = dg_dgrad_arho
    der[self.name][self.grad_arhou] = dg_dgrad_arhou
