from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters

class InterfaceClosuresParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerParameter("factory", "Factory")

class InterfaceClosures(object, metaclass=ABCMeta):
  def __init__(self, params):
    self.factory = params.get("factory")

  @abstractmethod
  def createAuxQuantities(self):
    pass
