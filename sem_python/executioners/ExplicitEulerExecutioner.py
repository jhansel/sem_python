from copy import deepcopy
import numpy as np

from .TransientExecutioner import TransientExecutioner, TransientExecutionerParameters

class ExplicitEulerExecutionerParameters(TransientExecutionerParameters):
  def __init__(self):
    TransientExecutionerParameters.__init__(self)

class ExplicitEulerExecutioner(TransientExecutioner):
  def __init__(self, params):
    TransientExecutioner.__init__(self, params)
    # create a mass matrix modified for strong constraints
    self.M_modified = deepcopy(self.M)
    self.applyStrongConstraintsToLinearSystemMatrix(self.M_modified)

    # invert mass matrix only once and keep it
    self.M_inv = np.linalg.inv(self.M_modified)

  def solve(self):
    # compute steady-state residual vector (as it would be on LHS)
    r_ss, J_ss = self.assembleSteadyStateSystemWithoutStrongConstraints(self.U_old, self.U_old)

    # compute linear system RHS vector
    b = np.matmul(self.M, self.U_old) - self.dt * r_ss

    # modify RHS for strong constraints
    self.applyStrongConstraintsToLinearSystemRHSVector(self.U_old, b)

    # compute the update
    self.U = np.matmul(self.M_inv, b)
