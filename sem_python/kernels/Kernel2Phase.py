from .Kernel import Kernel, KernelParameters

class Kernel2PhaseParameters(KernelParameters):
  def __init__(self):
    KernelParameters.__init__(self)
    self.registerIntParameter("phase", "Phase index (0 or 1)")
    self.registerParameter("var_enum", "Variable enumeration")

class Kernel2Phase(Kernel):
  def __init__(self, params):
    self.phase = params.get("phase")
    var_enum = params.get("var_enum")
    dof_handler = params.get("dof_handler")
    params.set("var_index", dof_handler.variable_index[var_enum][self.phase])
    Kernel.__init__(self, params)
    if self.phase == 0:
      self.sign = 1.0
    else:
      self.sign = -1.0

    # create list of relevant variable indices
    self.aA1_index = self.dof_handler.aA1_index[0]
    self.arhoA1_index = self.dof_handler.arhoA_index[0]
    self.arhouA1_index = self.dof_handler.arhouA_index[0]
    self.arhoEA1_index = self.dof_handler.arhoEA_index[0]
    self.arhoA2_index = self.dof_handler.arhoA_index[1]
    self.arhouA2_index = self.dof_handler.arhouA_index[1]
    self.arhoEA2_index = self.dof_handler.arhoEA_index[1]
    self.var_indices = [self.aA1_index, self.arhoA1_index, self.arhouA1_index, self.arhoEA1_index,
      self.arhoA2_index, self.arhouA2_index, self.arhoEA2_index]

    # create variable names
    self.vf1 = "vf1"
    self.arhoA1 = "arhoA1"
    self.arhouA1 = "arhouA1"
    self.arhoEA1 = "arhoEA1"
    self.arhoA2 = "arhoA2"
    self.arhouA2 = "arhouA2"
    self.arhoEA2 = "arhoEA2"

    self.rho1 = "rho1"
    self.u1 = "u1"
    self.E1 = "E1"
    self.v1 = "v1"
    self.e1 = "e1"
    self.p1 = "p1"
    self.T1 = "T1"

    self.rho2 = "rho2"
    self.u2 = "u2"
    self.E2 = "E2"
    self.v2 = "v2"
    self.e2 = "e2"
    self.p2 = "p2"
    self.T2 = "T2"
