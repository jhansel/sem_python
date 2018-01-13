import numpy as np

from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters


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
        f = data[self.vf] * data[self.visccoef_arhoA] * data["A"] * data[self.grad_rho] \
          + data[self.rho] * data[self.viscflux_aA]

        df_daA1 = (der[self.vf]["aA1"] * data[self.visccoef_arhoA] \
          + data[self.vf] * der[self.visccoef_arhoA]["aA1"]) * data["A"] * data[self.grad_rho] \
          + der[self.rho]["aA1"] * data[self.viscflux_aA] + data[self.rho] * der[self.viscflux_aA]["aA1"]
        df_darhoA1 = data[self.vf] * der[self.visccoef_arhoA]["arhoA1"] * data["A"] * data[self.grad_rho] \
          + der[self.rho]["arhoA1"] * data[self.viscflux_aA] + data[self.rho] * der[self.viscflux_aA]["arhoA1"]
        df_darhouA1 = data[self.vf] * der[self.visccoef_arhoA]["arhouA1"] * data["A"] * data[self.grad_rho] \
          + data[self.rho] * der[self.viscflux_aA]["arhouA1"]
        df_darhoEA1 = data[self.vf] * der[self.visccoef_arhoA]["arhoEA1"] * data["A"] * data[self.grad_rho] \
          + data[self.rho] * der[self.viscflux_aA]["arhoEA1"]
        df_darhoA2 = data[self.vf] * der[self.visccoef_arhoA]["arhoA2"] * data["A"] * data[self.grad_rho] \
          + der[self.rho]["arhoA2"] * data[self.viscflux_aA] + data[self.rho] * der[self.viscflux_aA]["arhoA2"]
        df_darhouA2 = data[self.vf] * der[self.visccoef_arhoA]["arhouA2"] * data["A"] * data[self.grad_rho] \
          + data[self.rho] * der[self.viscflux_aA]["arhouA2"]
        df_darhoEA2 = data[self.vf] * der[self.visccoef_arhoA]["arhoEA2"] * data["A"] * data[self.grad_rho] \
          + data[self.rho] * der[self.viscflux_aA]["arhoEA2"]
        df_dgrad_aA1 = data[self.vf] * data[self.visccoef_arhoA] * data["A"] * der[self.grad_rho]["grad_aA1"] \
          + data[self.rho] * der[self.viscflux_aA]["grad_aA1"]
        df_dgrad_arhoA = data[self.vf] * data[self.visccoef_arhoA] * data["A"] * der[self.grad_rho][
            self.grad_arhoA]

        data[self.name] = f
        der[self.name]["aA1"] = df_daA1
        der[self.name]["arhoA1"] = df_darhoA1
        der[self.name]["arhouA1"] = df_darhouA1
        der[self.name]["arhoEA1"] = df_darhoEA1
        der[self.name]["arhoA2"] = df_darhoA2
        der[self.name]["arhouA2"] = df_darhouA2
        der[self.name]["arhoEA2"] = df_darhoEA2
        der[self.name]["grad_aA1"] = df_dgrad_aA1
        der[self.name][self.grad_arhoA] = df_dgrad_arhoA
