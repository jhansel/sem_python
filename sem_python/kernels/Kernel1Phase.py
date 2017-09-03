from ..base.enums import ModelType
from .Kernel import Kernel, KernelParameters

class Kernel1PhaseParameters(KernelParameters):
  def __init__(self):
    KernelParameters.__init__(self)
    self.registerIntParameter("phase", "Phase index (0 or 1)")
    self.registerParameter("var_enum", "Variable enumeration")

class Kernel1Phase(Kernel):
  def __init__(self, params):
    self.phase = params.get("phase")
    self.var_enum = params.get("var_enum")
    dof_handler = params.get("dof_handler")
    params.set("var_index", dof_handler.variable_index[self.var_enum][self.phase])
    Kernel.__init__(self, params)

    # create list of relevant variable indices
    self.arho_index = self.dof_handler.arho_index[self.phase]
    self.arhou_index = self.dof_handler.arhou_index[self.phase]
    self.arhoE_index = self.dof_handler.arhoE_index[self.phase]
    if self.dof_handler.model_type == ModelType.TwoPhase:
      self.vf1_index = self.dof_handler.vf1_index[0]
      self.var_indices = [self.vf1_index, self.arho_index, self.arhou_index, self.arhoE_index]
    else:
      self.vf1_index = float("NaN")
      self.var_indices = [self.arho_index, self.arhou_index, self.arhoE_index]

    # create variable names
    phase_str = str(self.phase + 1)
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

    self.grad_arho = "grad_" + self.arho
    self.grad_arhou = "grad_" + self.arhou
    self.grad_arhoE = "grad_" + self.arhoE
