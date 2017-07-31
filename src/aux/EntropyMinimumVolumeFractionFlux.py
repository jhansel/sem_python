import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

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

  def compute(self, data, der):
    l = data[self.visccoef_vf] * data[self.grad_vf]
    dl_dlf_visc = data[self.grad_vf]
    dl_dgrad_vf = data[self.visccoef_vf]

    dl_dvf1 = dl_dlf_visc * der[self.visccoef_vf]["vf1"]
    dl_darho = dl_dlf_visc * der[self.visccoef_vf][self.arho]
    dl_darhou = dl_dlf_visc * der[self.visccoef_vf][self.arhou]
    dl_darhoE = dl_dlf_visc * der[self.visccoef_vf][self.arhoE]
    dl_dgrad_vf1 = dl_dgrad_vf * self.dgrad_vf_dgrad_vf1

    data[self.viscflux_vf] = l
    der[self.viscflux_vf] = {"vf1" : dl_dvf1, self.arho : dl_darho, self.arhou : dl_darhou,
      self.arhoE : dl_darhoE, "grad_vf1" : dl_dgrad_vf1}
