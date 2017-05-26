import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "src/kernels")
from Kernel import Kernel, KernelParameters

class Kernel1PhaseParameters(KernelParameters):
  def __init__(self):
    KernelParameters.__init__(self)
    self.registerIntParameter("phase", "Phase index (0 or 1)")

class Kernel1Phase(Kernel):
  def __init__(self, params, dof_handler, var_name):
    phase = params.get("phase")
    params.set("var_index", dof_handler.variable_index[var_name][phase])
    Kernel.__init__(self, params, dof_handler)

    # create list of relevant variable indices
    self.arho_index = self.dof_handler.arho_index[phase]
    self.arhou_index = self.dof_handler.arhou_index[phase]
    self.arhoE_index = self.dof_handler.arhoE_index[phase]
    if self.dof_handler.model_type == ModelType.TwoPhase:
      self.vf1_index = self.dof_handler.vf1_index[0]
      self.var_indices = [self.vf1_index, self.arho_index, self.arhou_index, self.arhoE_index]
    else:
      self.vf1_index = float("NaN")
      self.var_indices = [self.arho_index, self.arhou_index, self.arhoE_index]

    # create variable names
    phase_str = str(phase + 1)
    self.vf = "vf" + phase_str
    self.arho = "arho" + phase_str
    self.arhou = "arhou" + phase_str
    self.arhoE = "arhoE" + phase_str
    self.rho = "rho" + phase_str
    self.u = "u" + phase_str
    self.E = "E" + phase_str
    self.v = "v" + phase_str
    self.e = "e" + phase_str
    self.p = "p" + phase_str
