import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/executioners")
from TransientExecutioner import TransientExecutioner, TransientExecutionerParameters

sys.path.append(base_dir + "src/solvers")
from NonlinearSolver import NonlinearSolver

class ImplicitEulerExecutionerParameters(TransientExecutionerParameters):
  def __init__(self):
    TransientExecutionerParameters.__init__(self)

class ImplicitEulerExecutioner(TransientExecutioner):
  def __init__(self, params, model, ics, bcs, eos_map, interface_closures, gravity, dof_handler, mesh, nonlinear_solver_params, stabilization, factory):
    TransientExecutioner.__init__(self, params, model, ics, bcs, eos_map, interface_closures, gravity, dof_handler, mesh, nonlinear_solver_params, stabilization, factory)
    self.nonlinear_solver = NonlinearSolver(
      self.nonlinear_solver_params,
      self.assembleSystem,
      self.dof_handler)

  def assembleSystem(self, U):
    r_tr, J_tr = self.assembleTransientSystem(U)
    r_ss, J_ss = self.assembleSteadyStateSystemWithoutStrongBC(U)
    r = r_tr + self.dt * r_ss
    J = J_tr + self.dt * J_ss
    self.applyStrongBCNonlinearSystem(U, r, J)
    return (r, J)
