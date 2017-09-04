from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters

class AuxQuantityParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)

class AuxQuantity(object):
  __metaclass__ = ABCMeta

  def __init__(self, params):
    pass

  @abstractmethod
  def compute(self, data, der):
    pass
