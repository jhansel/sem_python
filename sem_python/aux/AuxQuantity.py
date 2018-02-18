from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class AuxQuantityParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)


class AuxQuantity(object, metaclass=ABCMeta):

    def __init__(self, params):
        pass

    @abstractmethod
    def compute(self, data, der):
        pass
