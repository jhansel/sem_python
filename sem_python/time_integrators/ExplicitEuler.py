from copy import deepcopy
import numpy as np

from .TimeIntegrator import TimeIntegrator, TimeIntegratorParameters


class ExplicitEulerParameters(TimeIntegratorParameters):

    def __init__(self, factory):
        TimeIntegratorParameters.__init__(self, factory)


class ExplicitEuler(TimeIntegrator):

    def __init__(self, params):
        TimeIntegrator.__init__(self, params)

        # create a mass matrix modified for strong constraints
        self.M_modified = deepcopy(self.M)
        self.assembly.applyConstraintsToLinearSystemMatrix(self.M_modified)

        # invert mass matrix only once and keep it
        self.M_inv = np.linalg.inv(self.M_modified)

    def performTimeStep(self, U):
        # compute steady-state residual vector (as it would be on LHS)
        r_ss, J_ss = self.assembly.assembleSteadyStateSystemWithoutConstraints(self.U_old)

        # compute linear system RHS vector
        b = np.matmul(self.M, self.U_old) - self.dt * r_ss

        # modify RHS for strong constraints
        self.assembly.applyConstraintsToLinearSystemRHSVector(self.U_old, b)

        # compute the update
        U = np.matmul(self.M_inv, b)
