from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class InterfaceClosuresParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("factory", "Factory")
        self.registerIntParameter("n_q", "Number of quadrature points per cell")


class InterfaceClosures(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.factory = params.get("factory")
        self.n_q = params.get("n_q")

    @abstractmethod
    def createAuxQuantities(self):
        pass
