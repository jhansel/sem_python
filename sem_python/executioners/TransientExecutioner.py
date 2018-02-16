from abc import abstractmethod
from copy import deepcopy
import numpy as np
from termcolor import colored

from ..base.enums import ModelType, VariableName
from .Executioner import Executioner, ExecutionerParameters
from ..utilities.error_utilities import error
from ..utilities.assembly_utilities import initializeDerivativeData


class TransientExecutionerParameters(ExecutionerParameters):

    def __init__(self, factory):
        ExecutionerParameters.__init__(self, factory)
        self.registerParameter("factory", "Factory")
        self.registerNamedSubblock("TimeStepSizer")
        self.registerBoolParameter("multiply_by_dt", "Multiply the nonlinear system by dt?", True)
        self.registerFloatParameter("ss_tol", "Tolerance for steady-state check")


class TransientExecutioner(Executioner):

    def __init__(self, params):
        Executioner.__init__(self, params)

        self.factory = params.get("factory")
        self.time_step_sizer = self.factory.createObjectOfType(params.get("TimeStepSizer"))
        self.multiply_by_dt = params.get("multiply_by_dt")

        if params.has("ss_tol"):
            self.check_ss = True
            self.ss_tol = params.get("ss_tol")
        else:
            self.check_ss = False

        self.U_old = deepcopy(self.U)

        # perform any setup particular to transient, such as assembling mass matrix
        self.assembly.performTransientSetup()

    @abstractmethod
    def solve(self):
        pass

    def run(self):
        while (self.time_step_sizer.transientIncomplete()):
            # compute time step size
            self.dt = self.time_step_sizer.getTimeStepSize(self.U)
            if self.verbose:
                self.time_step_sizer.printTimeStepInfo()

            # solve the time step
            self.solve()

            # check for steady-state
            if self.check_ss:
                dU_dt = (self.U - self.U_old) / self.dt
                dU_dt_norm = np.linalg.norm(dU_dt, 2)
                U_old_norm = np.linalg.norm(self.U_old, 2)
                U_change_norm = dU_dt_norm / U_old_norm
                if self.verbose:
                    print("Relative solution change: %e" % (U_change_norm))
                if U_change_norm < self.ss_tol:
                    if self.verbose:
                        print(colored("\nConverged to steady-state!\n", "green"))
                    return self.U

            # save old solution and increment time step index
            self.U_old = deepcopy(self.U)

        if self.verbose:
            print("")

        return self.U
