import numpy as np

from .TransientExecutioner import TransientExecutioner, TransientExecutionerParameters


class ImplicitEulerExecutionerParameters(TransientExecutionerParameters):

    def __init__(self, factory):
        TransientExecutionerParameters.__init__(self, factory)
        self.registerParameter("nonlinear_solver", "Nonlinear solver")


class ImplicitEulerExecutioner(TransientExecutioner):

    def __init__(self, params):
        TransientExecutioner.__init__(self, params)

        self.nonlinear_solver = params.get("nonlinear_solver")

    def assembleSystem(self, U):
        r_tr, J_tr = self.assembly.assembleTransientSystem(U, self.U_old)
        r_ss, J_ss = self.assembly.assembleSteadyStateSystemWithoutConstraints(U)
        if self.multiply_by_dt:
            r = r_tr + self.dt * r_ss
            J = J_tr + self.dt * J_ss
        else:
            r = r_tr / self.dt + r_ss
            J = J_tr / self.dt + J_ss
        self.assembly.applyConstraintsToNonlinearSystem(U, r, J)
        return (r, J)

    def solve(self):
        self.nonlinear_solver.solve(self.assembleSystem, self.U)
