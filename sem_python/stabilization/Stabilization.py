from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class StabilizationParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("factory", "Factory")
        self.registerParameter("dof_handler", "Degree of freedom handler")
        self.registerParameter("model_type", "Model type")
        self.registerParameter("quadrature", "Quadrature")


## Abstract base class for stabilization
class Stabilization(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.factory = params.get("factory")
        self.dof_handler = params.get("dof_handler")
        self.model_type = params.get("model_type")
        self.n_q = params.get("quadrature").n_q

    @abstractmethod
    def needSolutionGradients(self):
        pass

    @abstractmethod
    def createAuxQuantities(self):
        pass

    @abstractmethod
    def createIndependentPhaseKernels(self, phase):
        pass

    @abstractmethod
    def createPhaseInteractionKernels(self):
        pass
