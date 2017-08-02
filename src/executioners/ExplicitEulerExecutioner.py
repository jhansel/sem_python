from copy import deepcopy

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
    # create a mass matrix modified for Dirichlet BC
    self.M_modified = deepcopy(self.M)
    self.applyStrongBCLinearSystemMatrix(self.M_modified)

    # invert mass matrix only once and keep it
    self.M_inv = np.linalg.inv(self.M_modified)

  def solve(self):
    # compute steady-state residual vector (as it would be on LHS)
    r_ss, J_ss = self.assembleSteadyStateSystemWithoutStrongBC(self.U_old)

    # compute linear system RHS vector
    b = np.matmul(self.M, self.U_old) - self.dt * r_ss

    # modify RHS for Dirichlet BC
    self.applyStrongBCLinearSystemRHSVector(self.U_old, b)

    # compute the update
    self.U = np.matmul(self.M_inv, b)
