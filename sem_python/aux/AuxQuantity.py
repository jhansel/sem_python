from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class AuxQuantityParameters(Parameters):

    def __init__(self):
        Parameters.__init__(self)
        self.registerIntParameter("size", "Number of data values computed", 1)


class AuxQuantity(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.size = params.get("size")

    @abstractmethod
    def compute(self, data, der):
        pass
