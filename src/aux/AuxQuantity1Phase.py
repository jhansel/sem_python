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
    AuxQuantity.__init__(self)
    phase = str(params.get("phase") + 1)
    self.vf = "vf" + phase
    self.arho = "arho" + phase
    self.arhou = "arhou" + phase
    self.arhoE = "arhoE" + phase
    self.rho = "rho" + phase
    self.u = "u" + phase
    self.E = "E" + phase
    self.e = "e" + phase
    self.v = "v" + phase
    self.p = "p" + phase
    self.T = "T" + phase
    self.c = "c" + phase
