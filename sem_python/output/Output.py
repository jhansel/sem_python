from abc import ABCMeta, abstractmethod

from Parameters import Parameters

class OutputParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerParameter("dof_handler", "Degree of freedom handler")

class Output(object):
  __metaclass__ = ABCMeta

  def __init__(self, params):
    self.dof_handler = params.get("dof_handler")

  @abstractmethod
  def run(self, data):
    pass
