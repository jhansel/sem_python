import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/executioners")
from Executioner import Executioner, ExecutionerParameters

sys.path.append(base_dir + "src/solvers")
from NonlinearSolver import NonlinearSolver

class SteadyStateExecutionerParameters(ExecutionerParameters):
  def __init__(self):
    ExecutionerParameters.__init__(self)

class SteadyStateExecutioner(Executioner):
  def __init__(self, params):
    Executioner.__init__(self, params)

  def run(self):
    self.nonlinear_solver_params["assemble_system_function"] = self.assembleSteadyStateSystem
    self.nonlinear_solver_params["dof_handler"] = self.dof_handler
    nonlinear_solver = self.factory.createObject("NonlinearSolver", self.nonlinear_solver_params)
    nonlinear_solver.solve(self.U)

    return self.U
