import numpy as np

from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class EntropyMinimumMassFluxParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

## Entropy minimum viscous flux for the mass equation
#
# For phase \f$k\f$, this flux is computed as
# \f[
#   f_k \equiv \alpha_k\nu_k\nabla\rho_k + \rho_k l_k ,
# \f]
# where \f$l_k\f$ is the viscous flux for the volume fraction equation.
#
class EntropyMinimumMassFlux(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.viscflux_arho

  def compute(self, data, der):
    f = data[self.vf] * data[self.visccoef_arho] * data[self.grad_rho] \
      + data[self.rho] * data[self.viscflux_vf]

    df_dvf1 = (der[self.vf]["vf1"] * data[self.visccoef_arho] \
      + data[self.vf] * der[self.visccoef_arho]["vf1"]) * data[self.grad_rho] \
      + der[self.rho]["vf1"] * data[self.viscflux_vf] + data[self.rho] * der[self.viscflux_vf]["vf1"]
    df_darho = data[self.vf] * der[self.visccoef_arho][self.arho] * data[self.grad_rho] \
      + der[self.rho][self.arho] * data[self.viscflux_vf] + data[self.rho] * der[self.viscflux_vf][self.arho]
    df_darhou = data[self.vf] * der[self.visccoef_arho][self.arhou] * data[self.grad_rho] \
      + data[self.rho] * der[self.viscflux_vf][self.arhou]
    df_darhoE = data[self.vf] * der[self.visccoef_arho][self.arhoE] * data[self.grad_rho] \
      + data[self.rho] * der[self.viscflux_vf][self.arhoE]
    df_dgrad_vf1 = data[self.vf] * data[self.visccoef_arho] * der[self.grad_rho]["grad_vf1"] \
      + data[self.rho] * der[self.viscflux_vf]["grad_vf1"]
    df_dgrad_arho = data[self.vf] * data[self.visccoef_arho] * der[self.grad_rho][self.grad_arho]

    data[self.name] = f
    der[self.name] = {"vf1" : df_dvf1, self.arho : df_darho, self.arhou : df_darhou,
      self.arhoE : df_darhoE, "grad_vf1" : df_dgrad_vf1, self.grad_arho : df_dgrad_arho,
      self.grad_arhou: 0, self.grad_arhoE: 0}
