from .Executioner import Executioner, ExecutionerParameters
from ..solvers.NonlinearSolver import NonlinearSolver


class SteadyStateExecutionerParameters(ExecutionerParameters):

    def __init__(self, factory):
        ExecutionerParameters.__init__(self, factory)


class SteadyStateExecutioner(Executioner):

    def __init__(self, params):
        Executioner.__init__(self, params)

    def assembleSystem(self, U):
        r, J = self.assembleSteadyStateSystemWithoutStrongConstraints(U)
        self.applyStrongConstraintsToNonlinearSystem(U, r, J)
        return (r, J)

    def run(self):
        self.nonlinear_solver.solve(self.assembleSystem, self.U)

        return self.U
