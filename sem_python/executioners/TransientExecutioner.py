from termcolor import colored

from .Executioner import Executioner, ExecutionerParameters


class TransientExecutionerParameters(ExecutionerParameters):

    def __init__(self, factory):
        ExecutionerParameters.__init__(self, factory)
        self.registerParameter("factory", "Factory")
        self.registerNamedSubblock("TimeStepSizer")
        self.registerNamedSubblock("TimeIntegrator")
        self.registerFloatParameter("ss_tol", "Tolerance for steady-state check")


class TransientExecutioner(Executioner):

    def __init__(self, params):
        Executioner.__init__(self, params)

        self.factory = params.get("factory")
        self.time_step_sizer = self.factory.createObjectOfType(params.get("TimeStepSizer"))
        self.time_integrator = self.factory.createObjectOfType(params.get("TimeIntegrator"))

        if params.has("ss_tol"):
            self.check_ss = True
            self.ss_tol = params.get("ss_tol")
        else:
            self.check_ss = False

        # perform any setup particular to transient, such as assembling mass matrix
        self.assembly.performTransientSetup()

    def run(self):
        U = self.computeInitialSolution()

        while (self.time_step_sizer.transientIncomplete()):
            # compute time step size
            self.dt = self.time_step_sizer.getTimeStepSize(U)
            if self.verbose:
                self.time_step_sizer.printTimeStepInfo()

            # solve the time step
            self.time_integrator.takeTimeStep(U, self.dt)

            # check for steady-state
            if self.check_ss:
                U_change_norm = self.time_integrator.computeSteadyStateNorm(U)
                if self.verbose:
                    print("Relative solution change: %e" % (U_change_norm))
                if U_change_norm < self.ss_tol:
                    if self.verbose:
                        print(colored("\nConverged to steady-state!\n", "green"))
                    return U

        if self.verbose:
            print("")

        return U
