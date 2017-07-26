import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/kernels")
from Kernel import Kernel, KernelParameters

class Kernel2PhaseParameters(KernelParameters):
  def __init__(self):
    KernelParameters.__init__(self)
    self.registerIntParameter("phase", "Phase index (0 or 1)")

class Kernel2Phase(Kernel):
  def __init__(self, params, dof_handler, var_name):
    self.phase = params.get("phase")
    params.set("var_index", dof_handler.variable_index[var_name][self.phase])
    Kernel.__init__(self, params, dof_handler)
    if self.phase == 0:
      self.sign = 1.0
    else:
      self.sign = -1.0

    # create list of relevant variable indices
    self.vf1_index = self.dof_handler.vf1_index[0]
    self.arho1_index = self.dof_handler.arho_index[0]
    self.arhou1_index = self.dof_handler.arhou_index[0]
    self.arhoE1_index = self.dof_handler.arhoE_index[0]
    self.arho2_index = self.dof_handler.arho_index[1]
    self.arhou2_index = self.dof_handler.arhou_index[1]
    self.arhoE2_index = self.dof_handler.arhoE_index[1]
    self.var_indices = [self.vf1_index, self.arho1_index, self.arhou1_index, self.arhoE1_index,
      self.arho2_index, self.arhou2_index, self.arhoE2_index]

    # create variable names
    self.vf1 = "vf1"
    self.arho1 = "arho1"
    self.arhou1 = "arhou1"
    self.arhoE1 = "arhoE1"
    self.arho2 = "arho2"
    self.arhou2 = "arhou2"
    self.arhoE2 = "arhoE2"

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
