import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
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

  def compute(self, data, der):
    g = data[self.vf] * data[self.visccoef_arhou] * data[self.rho] * data[self.grad_u] \
      + data[self.u] * data[self.viscflux_arho]

    dg_dvf1 = (der[self.vf]["vf1"] * data[self.visccoef_arhou] * data[self.rho] \
      + data[self.vf] * der[self.visccoef_arhou]["vf1"] * data[self.rho] \
      + data[self.vf] * data[self.visccoef_arhou] * der[self.rho]["vf1"]) * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arho]["vf1"]
    dg_darho = data[self.vf] * (der[self.visccoef_arhou][self.arho] * data[self.rho] \
      + data[self.visccoef_arhou] * der[self.rho][self.arho]) * data[self.grad_u] \
      + der[self.u][self.arho] * data[self.viscflux_arho] + data[self.u] * der[self.viscflux_arho][self.arho]
    dg_darhou = data[self.vf] * der[self.visccoef_arhou][self.arhou] * data[self.rho] * data[self.grad_u] \
      + der[self.u][self.arhou] * data[self.viscflux_arho] + data[self.u] * der[self.viscflux_arho][self.arhou]
    dg_darhoE = data[self.vf] * der[self.visccoef_arhou][self.arhoE] * data[self.rho] * data[self.grad_u] \
      + data[self.u] * der[self.viscflux_arho][self.arhoE]
    dg_dgrad_vf1 = data[self.u] * der[self.viscflux_arho]["grad_vf1"]
    dg_dgrad_arho = data[self.vf] * data[self.visccoef_arhou] * data[self.rho] * der[self.grad_u][self.grad_arho] \
      + data[self.u] * der[self.viscflux_arho][self.grad_arho]
    dg_dgrad_arhou = data[self.vf] * data[self.visccoef_arhou] * data[self.rho] * der[self.grad_u][self.grad_arhou]

    data[self.viscflux_arhou] = g
    der[self.viscflux_arhou] = {"vf1" : dg_dvf1, self.arho : dg_darho, self.arhou : dg_darhou,
      self.arhoE : dg_darhoE, "grad_vf1" : dg_dgrad_vf1, self.grad_arho : dg_dgrad_arho,
      self.grad_arhou : dg_dgrad_arhou}
