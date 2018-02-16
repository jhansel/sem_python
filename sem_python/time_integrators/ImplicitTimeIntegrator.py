from abc import ABCMeta, abstractmethod

from .TimeIntegrator import TimeIntegrator, TimeIntegratorParameters


class ImplicitTimeIntegratorParameters(TimeIntegratorParameters):

    def __init__(self, factory):
        TimeIntegratorParameters.__init__(self, factory)
        self.registerParameter("nonlinear_solver", "Nonlinear solver")


class ImplicitTimeIntegrator(TimeIntegrator):

    def __init__(self, params):
        TimeIntegrator.__init__(self, params)
        self.nonlinear_solver = params.get("nonlinear_solver")

    ##
    # Returns a 2-tuple of the nonlinear residual and its Jacobian matrix
    #
    # @param[in] U   solution vector iterate
    # @returns 2-tuple of the nonlinear residual and its Jacobian matrix
    #
    @abstractmethod
    def assembleNonlinearSystem(self, U):
        pass

    def performTimeStep(self, U):
        self.nonlinear_solver.solve(self.assembleNonlinearSystem, U)
