from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class AssemblyParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)


class Assembly(object, metaclass=ABCMeta):

    def __init__(self, params):
        pass

    def performTransientSetup(self):
        pass

    @abstractmethod
    def assembleSteadyStateSystemWithoutConstraints(self, U):
        pass

    @abstractmethod
    def assembleTransientSystem(self, U):
        pass

    @abstractmethod
    def applyConstraintsToNonlinearSystem(self, U, r, J):
        pass

    @abstractmethod
    def applyConstraintsToLinearSystemMatrix(self, A):
        pass

    @abstractmethod
    def applyConstraintsToLinearSystemRHSVector(self, U, b):
        pass
