from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class TimeStepSizerParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerFloatParameter("end_time", "End time")


class TimeStepSizer(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.end_time = params.get("end_time")

        # tolerance to prevent small final time steps due to floating point precision error
        self.end_tolerance = 1e-12

        # current time
        self.t = 0.0

        # current time step size
        self.dt = None

        # time index
        self.t_index = 0

        # flag that transient is still incomplete
        self.transient_incomplete = True

    ##
    # Returns whether the transient is still incomplete
    #
    # @returns whether the transient is still incomplete
    #
    def transientIncomplete(self):
        return self.transient_incomplete

    ##
    # Prints time step information
    #
    def printTimeStepInfo(self):
        print("\nTime step %i: t = %g, dt = %g" % (self.t_index, self.t, self.dt))

    ##
    # Returns the time step size and updates transient information
    #
    # @param[in] U   solution vector
    # @returns time step size
    #
    def getTimeStepSize(self, U):
        # get the time step size from the derived class
        self.dt = self.getTimeStepSizeInternal(U)

        # decrease time step size if at end of transient
        if (self.t + self.dt + self.end_tolerance >= self.end_time):
            self.transient_incomplete = False
            self.dt = self.end_time - self.t

        # update time and time index
        self.t += self.dt
        self.t_index += 1

        return self.dt

    ##
    # Returns the time step size
    #
    # @param[in] U   solution vector
    # @returns time step size
    #
    @abstractmethod
    def getTimeStepSizeInternal(self, U):
        pass
