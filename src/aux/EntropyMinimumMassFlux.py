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
    self.name = self.viscflux_arhoA

  def compute(self, data, der):
    f = data[self.vf] * data[self.visccoef_arhoA] * data[self.grad_rho] \
      + data[self.rho] * data[self.viscflux_vf]

    df_dvf1 = (der[self.vf]["vf1"] * data[self.visccoef_arhoA] \
      + data[self.vf] * der[self.visccoef_arhoA]["vf1"]) * data[self.grad_rho] \
      + der[self.rho]["vf1"] * data[self.viscflux_vf] + data[self.rho] * der[self.viscflux_vf]["vf1"]
    df_darhoA1 = data[self.vf] * der[self.visccoef_arhoA]["arhoA1"] * data[self.grad_rho] \
      + der[self.rho]["arhoA1"] * data[self.viscflux_vf] + data[self.rho] * der[self.viscflux_vf]["arhoA1"]
    df_darhouA1 = data[self.vf] * der[self.visccoef_arhoA]["arhouA1"] * data[self.grad_rho] \
      + data[self.rho] * der[self.viscflux_vf]["arhouA1"]
    df_darhoEA1 = data[self.vf] * der[self.visccoef_arhoA]["arhoEA1"] * data[self.grad_rho] \
      + data[self.rho] * der[self.viscflux_vf]["arhoEA1"]
    df_darhoA2 = data[self.vf] * der[self.visccoef_arhoA]["arhoA2"] * data[self.grad_rho] \
      + der[self.rho]["arhoA2"] * data[self.viscflux_vf] + data[self.rho] * der[self.viscflux_vf]["arhoA2"]
    df_darhouA2 = data[self.vf] * der[self.visccoef_arhoA]["arhouA2"] * data[self.grad_rho] \
      + data[self.rho] * der[self.viscflux_vf]["arhouA2"]
    df_darhoEA2 = data[self.vf] * der[self.visccoef_arhoA]["arhoEA2"] * data[self.grad_rho] \
      + data[self.rho] * der[self.viscflux_vf]["arhoEA2"]
    df_dgrad_vf1 = data[self.vf] * data[self.visccoef_arhoA] * der[self.grad_rho]["grad_vf1"] \
      + data[self.rho] * der[self.viscflux_vf]["grad_vf1"]
    df_dgrad_arhoA = data[self.vf] * data[self.visccoef_arhoA] * der[self.grad_rho][self.grad_arhoA]

    data[self.name] = f
    der[self.name]["vf1"] = df_dvf1
    der[self.name]["arhoA1"] = df_darhoA1
    der[self.name]["arhouA1"] = df_darhouA1
    der[self.name]["arhoEA1"] = df_darhoEA1
    der[self.name]["arhoA2"] = df_darhoA2
    der[self.name]["arhouA2"] = df_darhouA2
    der[self.name]["arhoEA2"] = df_darhoEA2
    der[self.name]["grad_vf1"] = df_dgrad_vf1
    der[self.name][self.grad_arhoA] = df_dgrad_arhoA
