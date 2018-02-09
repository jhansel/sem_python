import numpy as np

from .TransientExecutioner import TransientExecutioner, TransientExecutionerParameters


class ImplicitEulerExecutionerParameters(TransientExecutionerParameters):

    def __init__(self, factory):
        TransientExecutionerParameters.__init__(self, factory)


class ImplicitEulerExecutioner(TransientExecutioner):

    def __init__(self, params):
        TransientExecutioner.__init__(self, params)

    def assembleSystem(self, U):
        r_tr, J_tr = self.assembleTransientSystem(U)
        r_ss, J_ss = self.assembleSteadyStateSystemWithoutConstraints(U)
        if self.multiply_by_dt:
            r = r_tr + self.dt * r_ss
            J = J_tr + self.dt * J_ss
        else:
            r = r_tr / self.dt + r_ss
            J = J_tr / self.dt + J_ss
        self.applyConstraintsToNonlinearSystem(U, r, J)
        return (r, J)

    def solve(self):
        self.nonlinear_solver.solve(self.assembleSystem, self.U)
