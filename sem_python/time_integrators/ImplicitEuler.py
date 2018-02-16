from .ImplicitTimeIntegrator import ImplicitTimeIntegrator, ImplicitTimeIntegratorParameters

class ImplicitEulerParameters(ImplicitTimeIntegratorParameters):

    def __init__(self, factory):
        ImplicitTimeIntegratorParameters.__init__(self, factory)


class ImplicitEuler(ImplicitTimeIntegrator):

    def __init__(self, params):
        ImplicitTimeIntegrator.__init__(self, params)

    def assembleNonlinearSystem(self, U):
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
