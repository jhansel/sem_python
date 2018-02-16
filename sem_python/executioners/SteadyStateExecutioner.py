from .Executioner import Executioner, ExecutionerParameters


class SteadyStateExecutionerParameters(ExecutionerParameters):

    def __init__(self, factory):
        ExecutionerParameters.__init__(self, factory)
        self.registerParameter("nonlinear_solver", "Nonlinear solver")


class SteadyStateExecutioner(Executioner):

    def __init__(self, params):
        Executioner.__init__(self, params)
        self.nonlinear_solver = params.get("nonlinear_solver")

    def assembleSystem(self, U):
        r, J = self.assembly.assembleSteadyStateSystemWithoutConstraints(U)
        self.assembly.applyConstraintsToNonlinearSystem(U, r, J)
        return (r, J)

    def run(self):
        U = self.computeInitialSolution()
        self.nonlinear_solver.solve(self.assembleSystem, U)

        return U
