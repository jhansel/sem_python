import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/executioners")
from TransientExecutioner import TransientExecutioner, TransientExecutionerParameters

class ExplicitEulerExecutionerParameters(TransientExecutionerParameters):
  def __init__(self):
    TransientExecutionerParameters.__init__(self)

class ExplicitEulerExecutioner(TransientExecutioner):
  def __init__(self, params, model_type, ics, bcs, eos_map, interface_closures, gravity, dof_handler, mesh, nonlinear_solver_params, stabilization, factory):
    TransientExecutioner.__init__(self, params, model_type, ics, bcs, eos_map, interface_closures, gravity, dof_handler, mesh, nonlinear_solver_params, stabilization, factory)

  def assembleSystem(self, U):
    r_tr, J_tr = self.assembleTransientSystem(U)
    r_ss, J_ss = self.assembleSteadyStateSystemWithoutStrongBC(self.U_old)
    r = r_tr + self.dt * r_ss
    J = J_tr
    self.applyStrongBC(U, r, J)
    return (r, J)
