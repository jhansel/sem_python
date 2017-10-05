from .Executioner import Executioner, ExecutionerParameters
from ..solvers.NonlinearSolver import NonlinearSolver

class SteadyStateExecutionerParameters(ExecutionerParameters):
  def __init__(self):
    ExecutionerParameters.__init__(self)

class SteadyStateExecutioner(Executioner):
  def __init__(self, params):
    Executioner.__init__(self, params)

  def assembleSystem(self, U):
    r, J = self.assembleSteadyStateSystemWithoutStrongConstraints(U)
    self.applyStrongConstraintsToNonlinearSystem(U, r, J)
    return (r, J)

  def run(self):
    self.nonlinear_solver_params["assemble_system_function"] = self.assembleSystem
    self.nonlinear_solver_params["dof_handler"] = self.dof_handler
    nonlinear_solver = self.factory.createObject("NonlinearSolver", self.nonlinear_solver_params)
    nonlinear_solver.solve(self.U)

    return self.U
