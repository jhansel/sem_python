import numpy as np

from .TransientExecutioner import TransientExecutioner, TransientExecutionerParameters

class ImplicitEulerExecutionerParameters(TransientExecutionerParameters):
  def __init__(self):
    TransientExecutionerParameters.__init__(self)

class ImplicitEulerExecutioner(TransientExecutioner):
  def __init__(self, params):
    TransientExecutioner.__init__(self, params)
    self.nonlinear_solver_params["assemble_system_function"] = self.assembleSystem
    self.nonlinear_solver_params["dof_handler"] = self.dof_handler
    self.nonlinear_solver = self.factory.createObject("NonlinearSolver", self.nonlinear_solver_params)

  def assembleSystem(self, U):
    r_tr, J_tr = self.assembleTransientSystem(U)
    r_ss, J_ss = self.assembleSteadyStateSystemWithoutStrongBC(U)
    r = r_tr + self.dt * r_ss
    J = J_tr + self.dt * J_ss
    self.applyStrongConstraintsToNonlinearSystem(U, r, J)
    return (r, J)
