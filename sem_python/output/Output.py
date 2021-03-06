from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters


class OutputParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("dof_handler", "Degree of freedom handler")


class Output(object, metaclass=ABCMeta):

    def __init__(self, params):
        self.dof_handler = params.get("dof_handler")

    @abstractmethod
    def run(self, data):
        pass
