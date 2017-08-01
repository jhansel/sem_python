import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity import AuxQuantity, AuxQuantityParameters

class AuxQuantity1PhaseParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)
    self.registerIntParameter("phase", "Phase index (0 or 1)")

class AuxQuantity1Phase(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self, params)
    self.phase_int = params.get("phase")
    self.phase = str(self.phase_int + 1)
    if self.phase_int == 0:
      self.sign = 1.0
      self.dgrad_vf_dgrad_vf1 = 1.0
    else:
      self.sign = -1.0
      self.dgrad_vf_dgrad_vf1 = -1.0

    self.vf = "vf" + self.phase
    self.arho = "arho" + self.phase
    self.arhou = "arhou" + self.phase
    self.arhoE = "arhoE" + self.phase
    self.rho = "rho" + self.phase
    self.rhoe = "rhoe" + self.phase
    self.u = "u" + self.phase
    self.E = "E" + self.phase
    self.e = "e" + self.phase
    self.v = "v" + self.phase
    self.p = "p" + self.phase
    self.T = "T" + self.phase
    self.c = "c" + self.phase

    self.grad_arho = "grad_" + self.arho
    self.grad_arhou = "grad_" + self.arhou
    self.grad_arhoE = "grad_" + self.arhoE
    self.grad_vf = "grad_" + self.vf
    self.grad_rho = "grad_" + self.rho
    self.grad_u = "grad_" + self.u
    self.grad_rhoe = "grad_" + self.rhoe

    self.visccoef_vf = "visccoef_" + self.vf
    self.visccoef_arho = "visccoef_" + self.arho
    self.visccoef_arhou = "visccoef_" + self.arhou
    self.visccoef_arhoE = "visccoef_" + self.arhoE

    self.viscflux_vf = "viscflux_" + self.vf
    self.viscflux_arho = "viscflux_" + self.arho
    self.viscflux_arhou = "viscflux_" + self.arhou
    self.viscflux_arhoE = "viscflux_" + self.arhoE
