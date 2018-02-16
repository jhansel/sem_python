from abc import ABCMeta, abstractmethod
from copy import deepcopy
import numpy as np

from ..input.Parameters import Parameters


class TimeIntegratorParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerBoolParameter("multiply_by_dt", "Multiply nonlinear system by time step size", True)
        self.registerParameter("assembly", "Assembly")


class TimeIntegrator(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.multiply_by_dt = params.get("multiply_by_dt")
        self.assembly = params.get("assembly")

        self.U_old = None
        self.dt_old = None

    def takeTimeStep(self, U, dt):
        self.U_old = deepcopy(U)
        self.dt = dt

        self.performTimeStep(U)

        self.dt_old = dt

    @abstractmethod
    def performTimeStep(self, U):
        pass

    def computeSteadyStateNorm(self, U):
        dU_dt = (U - self.U_old) / self.dt_old
        dU_dt_norm = np.linalg.norm(dU_dt, 2)
        U_old_norm = np.linalg.norm(self.U_old, 2)

        return dU_dt_norm / U_old_norm
